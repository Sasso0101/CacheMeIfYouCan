#!/bin/bash
set -e

rm -rf build
cmake -DDBG_FRONTIER_SIZE=OFF -DDBG_THREAD_BALANCE=OFF -B build
cd build
make BFS # BFS_REF
cd ..
# sleep 2

export OMP_NUM_THREADS=24
export TIMEOUT=60
export NITER=20

# <num_runs> <timeout> <log_filename> <executable> <executable_args>
# scripts/run_on_all_datasets.sh 10 $TIMEOUT runtime_sequential_naive build/BFS_REF

for threads in 1 6 12 24 32 48 64 96; do
    export OMP_NUM_THREADS=$threads
    scripts/run_on_all_datasets.sh $NITER $TIMEOUT "runtime_classic$OMP_NUM_THREADS" build/BFS classic
    scripts/run_on_all_datasets.sh $NITER $TIMEOUT "runtime_small$OMP_NUM_THREADS" build/BFS small
    scripts/run_on_all_datasets.sh $NITER $TIMEOUT "runtime_large$OMP_NUM_THREADS" build/BFS large
done


# scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_classic$OMP_NUM_THREADS" build/BFS classic
# scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_small$OMP_NUM_THREADS" build/BFS small
# scripts/run_on_all_datasets.sh 10 $TIMEOUT "runtime_large$OMP_NUM_THREADS" build/BFS large