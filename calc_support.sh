#!/bin/bash

# check BS convergence and compute support
# NOTE: this script must be called from "<LONGREADS_ROOT>/trees/ref_<TREE>.fasta/bootstrap" dir

raxml=/home/alexey/hits/raxml-git/raxmlHPC-AVX2

cat ./*/bootstrap.raxml.bootstraps > bsrep.tre

rm RAxML_*.autoMRE

#check convergence
$raxml -n autoMRE -z bsrep.tre -I autoMRE -m GTRGAMMA -p 3 -B 0.03 --bootstop-perms=1000

rm RAxML_*.support

# map BS support values on the ML tree
$raxml -n support -z bsrep.tre -t ../best.tre -m GTRGAMMA -p 3 -f b

rm bsrep.tre
