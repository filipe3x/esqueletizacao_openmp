#!/bin/sh
#
#PBS -N ske_comm
#PBS -l walltime=04:59:59
#PBS -l nodes=8:r641:ppn=32:myri

#PBS -m abe
#PBS -M luis.fonseca156@gmail.com
#PBS -M filipe3x@hotmail.com

#PBS -j oe

MYWORKPLACE="/home/$USER/PCP/esqueletizacao_openmp/"

job="jobscript.sh"

log="results$(date +%d-%m-%Y-%H%M%S).log"

cd $MYWORKPLACE

#source ./load_openmp.sh
source ./load_mpi.sh

mkdir -p results

./$job > results/$log

