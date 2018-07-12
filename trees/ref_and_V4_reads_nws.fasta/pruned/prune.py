#!/usr/bin/python

from ete3 import Tree

BASE="/home/lucas/Projects/long-reads"

IN=BASE + "/trees/ref_and_V4_reads_nws.fasta/best.tre"
NAMES=BASE + "/data/ref_names.txt"
with open(NAMES) as f:
    TIPS = f.read().splitlines()
OUT=BASE + "/trees/ref_and_V4_reads_nws.fasta/pruned/pruned.newick"

t = Tree( IN )
t.prune(TIPS)
t.unroot()
t.write(format=5, outfile=OUT)
