#!/bin/bash
#@ job_name = longreads_search
#@ job_type = parallel
#@ energy_policy_tag = raxng_energy_tag
#@ minimize_time_to_solution = yes
#@ class = micro
#@ node = 1
#@ island_count = 1,1
#@ tasks_per_node = 8
#@ network.MPI = sn_all,not_shared,us
#@ output = job.$(schedd_host).$(jobid).out
#@ error =  job.$(schedd_host).$(jobid).err
#@ notification=always
#@ notify_user=pierre.barbera@h-its.org
#@ wall_clock_limit = 48:00:00
#@ queue

#SBATCH -o job.long_reads_search.%j.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -B 2:8:1
#SBATCH --threads-per-core=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00

source ~/haswell_modules.sh

THREADS=40

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ALI=${BASE}/data/ref_only_nws.fasta
TDIR=${BASE}/trees/ref_only_nws.fasta/UNCONS_RAND
TREE=${TDIR}/search.raxml.bestTree
MODELFILE=${TDIR}/search.raxml.bestModel
QRY=$1

if [ -d ${BASE} -a -r ${ALI} -a -r ${QRY} -a -r ${TREE} -a -r ${MODELFILE} ]
then echo "all paths OK!"
else echo "ERROR: a path was wrong!" && exit
fi

# extract model string
MODEL="$(head -n 1 ${MODELFILE} | awk -F "," '{print $1}')"

QRYNAME=$(basename ${QRY})

OUT=${BASE}/jplace/${QRYNAME}/

mkdir -p ${OUT}
rm -f ${OUT}/*

echo "start at `date`"

epa-ng --tree ${TREE} --ref-msa ${ALI} --query ${QRY} --outdir ${OUT} --model ${MODEL} --threads ${THREADS} --no-heur --no-pre-mask --verbose

echo "end at `date`"



