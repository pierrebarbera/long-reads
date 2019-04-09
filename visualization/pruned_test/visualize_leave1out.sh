#!/bin/bash

BASE=/home/lucas/Projects/long-reads
GENESIS=/home/lucas/Dropbox/HITS/genesis/bin/apps

mkdir ${BASE}/visualization/pruned/long
mkdir ${BASE}/visualization/pruned/V4

${GENESIS}/long_read_eval ${BASE}/jplace/pruned_long/long_read_queries/epa_result.jplace ${BASE}/visualization/pruned/long
${GENESIS}/long_read_eval ${BASE}/jplace/pruned_V4/V4_read_queries/epa_result.jplace     ${BASE}/visualization/pruned/V4
