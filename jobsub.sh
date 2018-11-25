#!/bin/sh
#
#PBS -N skeletonize
#PBS -l walltime=23:59:59
#PBS -l nodes=1:r641:ppn=32

#PBS -m abe
#PBS -M filipe3x@hotmail.com

module load papi/5.4.1

export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=false
export OMP_PROC_BIND=true

MYWORKPLACE="/home/a57812/PCP/esqueletizacao_openmp/"

job="jobscript.sh"

log="results$(date +%d-%m-%Y-%H%M%S).log"

cd $MYWORKPLACE

mkdir -p results

./$job > results/$log



