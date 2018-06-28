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

BASE=/hits/basement/sco/barbera/long-reads
ALI=${1:-${BASE}/data/ref_and_V4_reads_nws.fasta}

NAME[0]=UNCONS

NUM=${3:-10}

STARTF="--tree pars{${NUM}}"; NAME[1]=PARS

#for mod in "$@"; do

if [ $2 = rand ];
then STARTF="--tree rand{${NUM}}"; NAME[1]=RAND
fi

#done

ALINAME=$(basename ${ALI})

NAME=$(printf "_%s" "${NAME[@]}")
NAME=${NAME:1}

OUT=${BASE}/trees/${ALINAME}/${NAME}

mkdir -p ${OUT}
# rm -f ${OUT}/*

if [ -d ${BASE} -a -d ${OUT} -a -r ${ALI} ]
then echo "all paths OK!"
else echo "ERROR: a path was wrong!" && exit
fi

echo "start at `date`"

raxml-ng --search --msa ${ALI} --model GTR+G --prefix ${OUT}/search ${STARTF} --threads 8

echo "end at `date`"



