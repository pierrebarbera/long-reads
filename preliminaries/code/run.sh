#!/bin/bash

basedir=long-reads
genesis=genesis/bin/apps

rm *.svg
rm *.nexus

${genesis}/long_read_eval ${basedir}/RAxML_leaveOneOutResults.l1out_seq_blo.jplace ${basedir}
