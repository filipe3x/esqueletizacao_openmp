#!/bin/sh
#
#PBS -N skeletonize
#PBS -l walltime=23:59:59
#PBS -l nodes=1:r641:ppn=32

#PBS -m abe
#PBS -M luis.fonseca156@gmail.com
#PBS -M filipe3x@hotmail.com

MYWORKPLACE="/home/a57812/PCP/esqueletizacao_openmp/"

job="jobscript.sh"

log="results$(date +%d-%m-%Y-%H%M%S).log"

cd $MYWORKPLACE

#source ./load_openmp.sh
source ./load_mpi.sh

mkdir -p results

./$job > results/$log

