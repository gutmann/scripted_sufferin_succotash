#!/bin/sh
if [[ $1 ]]; then
    nnodes=$1
else
    nnodes=1
fi

if [[ $2 ]]; then
    hours=$2
else
    hours=4
fi

let time=${hours}*3600

echo "qsub -I -l select=${nnodes}:mpiprocs=36:ncpus=36:ompthreads=36 -l walltime=${time} -A P48500028 -q regular"
qsub -I -l select=${nnodes}:mpiprocs=36:ncpus=36:ompthreads=36 -l walltime=${time} -A P48500028 -q regular
