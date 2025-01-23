#!/bin/bash
set -e

cmake -DDBG_FRONTIER_SIZE=OFF -DDBG_THREAD_BALANCE=OFF -B build
cd build
make BFS BFS_REF
cd ..
sleep 2

export OMP_NUM_THREADS=24
export TIMEOUT=60

# <num_runs> <timeout> <log_filename> <executable> <executable_args>
scripts/run_on_all_datasets.sh 10 $TIMEOUT runtime_reference build/BFS_REF

scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_auto$OMP_NUM_THREADS" build/BFS

scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_small$OMP_NUM_THREADS" build/BFS small
scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_large$OMP_NUM_THREADS" build/BFS large
scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_verylarge$OMP_NUM_THREADS" build/BFS very_large