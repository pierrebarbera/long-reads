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

THREADS=10

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TAXONS=${BASE}/data/taxonomy.txt

for treename in "V4" "long" "ref_only"; do

  for queryname in "V4" "long"; do

    if [[ "${treename}" == "ref_only" ]];
      then  TPREFIX=""
      else  TPREFIX="pruned_"
    fi

    QRY=${BASE}/jplace/${TPREFIX}${treename}/${queryname}_read_queries/epa_result.jplace

    if [ -d ${BASE} -a -r ${QRY} -a -r ${TAXONS} ]
    then echo "all paths OK!"
    else echo "ERROR: a path was wrong!" && exit
    fi

    OUT=${BASE}/assign/${TPREFIX}${treename}/${queryname}_read_queries/

    mkdir -p ${OUT}
    rm -f ${OUT}/*

    gappa analyze assign --threads ${THREADS} --jplace-path ${QRY} --krona --taxon-file ${TAXONS} --out-dir ${OUT}

  done

done

echo "end at `date`"
