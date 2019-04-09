#!/usr/bin/python

# libraries
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
import seaborn as sns
import os

#######################################################################
#    Plot Functions
#######################################################################

# Plot the exact curve instead of a histogram approximation
def plot_exact_hist( data ):
	# Load list of distances, sort them
	data = np.sort( data )

	# For each distance, add a step on the y-axis,
	# normalized by how many values there are.
	steps = np.array(range(len(data))) / float(len(data))

	# We want to avoid the ugly jump in the beginning.
	# Set all steps in the beginning to the height at the last 0 step.
	last_zero_idx=0
	for i in range(len(data)):
		if data[i] == 0:
			last_zero_idx=i
		else:
			break
	for i in range(last_zero_idx):
		steps[i] = steps[last_zero_idx]

	# Plot the data
	plt.plot(data, steps)

def plot_dataset(datasets, measure):

	# Plot the graph.
	for dataset in datasets:
		filename = dataset + "/" + measure + "_list.csv"
		table = pd.read_csv( filename )
		
		# Get max x
		if measure == "edpl":
			max_x = 0.1
			data = table["EDPL"]
		elif measure == "lwr":
			max_x = 1.0
			data = table["LWR 1"]
		else:
			raise error
		plot_exact_hist( data )

	# Set texts.
	plt.legend(datasets, loc='lower right')
	plt.xlabel(measure)
	plt.ylabel("Cumulative Frequency")

	# Set 10 ticks on x and y axis, and set the desired axis ranges.
	plt.locator_params( nbins=10, axis='x' )
	plt.xlim( 0.0, max_x )
	plt.locator_params( nbins=5, axis='y' )
	# plt.ylim( 0.5, 1.0 )

	# Set y tick labels every 20% and grid lines every 10%.
	# Rename labels to be percent.
	ax = plt.gca()
	# ax.set_yticks( np.linspace( 0.0, 1.0, 6 ), minor=False)
	# ax.set_yticks( np.linspace( 0.1, 0.9, 5 ), minor=True)
	# ax.set_yticks( np.linspace( 0.6, 1.0, 3 ), minor=False)
	# ax.set_yticks( np.linspace( 0.7, 0.9, 2 ), minor=True)

	ax.yaxis.grid(True, which='major')
	ax.yaxis.grid(True, which='minor')
	ax.set_yticklabels([ '{:3.0f}%'.format(x*100) for x in ax.get_yticks() ])
	#ax.set_xticklabels([ '< {:1.0f}'.format(x) for x in ax.get_xticks() ], rotation=45)

	plt.tight_layout()
	#plt.figure().tight_layout()

	# Hide major tick labels and make them appear at the mid of the bins for histograms
	# if method == "histogram":
	# 	ax.set_xticklabels('')
	# 	ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5], minor=True)
	# 	ax.set_xticklabels(['0','1','2','3','4','5','6','7','8','9'], minor=True)
	# 	ax.xaxis.grid(False, which='minor')

	# Save files
	plt.savefig( "figures_png/" + measure + ".png", format='png' )
	plt.savefig( "figures_svg/" + measure + ".svg", format='svg' )
	#plt.show()

#######################################################################
#    Dataset Information
#######################################################################

datasets    = [ "long", "V4" ]
measure     = [ "edpl", "lwr" ]

#######################################################################
#    Plotting Loop
#######################################################################

#sns.set(font_scale=1.6)
plotnum = 1

# if not os.path.exists("figures_pdf"):
#     os.makedirs("figures_pdf")
if not os.path.exists("figures_png"):
    os.makedirs("figures_png")
if not os.path.exists("figures_svg"):
    os.makedirs("figures_svg")

# Make all needed combinations of plots
for measure in measure:
	print "Figure", plotnum, ":", measure
	plt.figure( plotnum )

	plot_dataset(datasets, measure)

	plt.close( plotnum )
	plotnum += 1

#plt.show()
