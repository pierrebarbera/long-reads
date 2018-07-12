#!/bin/bash

BASE=/home/lucas/Projects/long-reads/trees

# Three trees
RT=${BASE}/ref_only_nws.fasta/pruned/pruned.newick
LT=${BASE}/ref_and_long_reads_nws.fasta/pruned/pruned.newick
VT=${BASE}/ref_and_V4_reads_nws.fasta/pruned/pruned.newick

# Concat them
cat ${RT}  > ${BASE}/rf-dists/trees.nw
cat ${LT} >> ${BASE}/rf-dists/trees.nw
cat ${VT} >> ${BASE}/rf-dists/trees.nw

# RF dist
#raxml -m GTRCAT -z ${BASE}/rf-dists/trees.nw -f r -n rfd

# Visualize differences
NC=/home/lucas/Dropbox/HITS/genesis/bin/apps/newick_compare
${NC} ${RT} ${VT}
