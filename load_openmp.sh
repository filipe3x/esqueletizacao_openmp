#!/bin/sh
#

export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=false
export OMP_PROC_BIND=true

module load gcc/5.3.0
module load papi/5.4.1
