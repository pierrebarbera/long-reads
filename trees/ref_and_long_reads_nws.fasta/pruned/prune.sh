#!/bin/bash

BASE=/home/lucas/Projects/long-reads
NT=/home/lucas/Software/newick-tools/src/newick-tools

IN=${BASE}/trees/ref_and_long_reads_nws.fasta/best.tre
NAMES=${BASE}/data/long_reads_names.txt
TIPS=`cat ${NAMES} | sed 's/^\(.*\)$/ --prune_tips \1 /g'`
OUT=${BASE}/trees/ref_and_long_reads_nws.fasta/pruned/pruned.newick

#echo $TIPS

${NT} --tree_file ${IN} --info
#${NT} --tree_file ${IN} --output_file ${OUT} ${TIPS}
${NT} --tree_file ${OUT} --info

