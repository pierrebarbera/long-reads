#!/bin/bash

READLEN=$1

# number of the first QUERY sequence -> where to start leave-one-out test
L1OUT_START=1662

REPO_HOME=/hits/basement/sco/kozlov/long-reads/repo
RAXML=/hits/basement/sco/kozlov/sativa/raxml/raxmlHPC8-AVX2.PTHREADS
ALI=$REPO_HOME/data/ref_and_${READLEN}_reads_nws_sorted.fasta
ORIG_TREE=$REPO_HOME/trees/ref_and_${READLEN}_reads_nws.fasta/best.tre

# optimize branches/model and save to the RAxML binary file
$RAXML -T 16 -s $ALI -t $ORIG_TREE --no-seq-check -f e -m GTRGAMMAX -n eval_${READLEN} --verbose -p 1

# perform the leave-one-out test
$RAXML -T 16 -s $ALI -t RAxML_result.eval_${READLEN} --no-seq-check -f O --epa-accumulated-threshold 0.99 -m GTRGAMMAX -n l1out_${READLEN} --verbose -p 1 --l1out-start $L1OUT_START -R RAxML_binaryModelParameters.eval_${READLEN}

# generate alignment with a dummy sequence to place (will be used at the next step)
cp $ALI epalbl.fa

echo '>QQQQQ' >> epalbl.fa

tail -n 1 $ALI >> epalbl.fa

# run dummy placement to obtain a properly labelled reference tree
$RAXML -T 16 -s epalbl.fa -t RAxML_result.eval_${READLEN} --no-seq-check -f y --epa-accumulated-threshold 1.0 -m GTRGAMMA -n epalbl_${READLEN} --verbose -p 1 

rm epalbl.fa

# add reference tree into jplace file
grep '"tree":' RAxML_portableTree.epalbl_${READLEN}.jplace > reftree.json

sed '1r reftree.json' RAxML_leaveOneOutResults.l1out_${READLEN}.jplace > L1OUT_RESULT_${READLEN}.jplace

rm reftree.json
