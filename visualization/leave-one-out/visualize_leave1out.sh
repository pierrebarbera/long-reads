#!/bin/bash

BASE=/home/lucas/Projects/long-reads
GENESIS=/home/lucas/Dropbox/HITS/genesis/bin/apps

mkdir ${BASE}/visualization/long
mkdir ${BASE}/visualization/V4

${GENESIS}/long_read_eval ${BASE}/leave1out/L1OUT_RESULT_long.jplace ${BASE}/visualization/long
${GENESIS}/long_read_eval ${BASE}/leave1out/L1OUT_RESULT_V4.jplace   ${BASE}/visualization/V4
