#!/bin/bash
set -e

rm -rf build
rm -f sources.txt
cmake -DDBG_FRONTIER_SIZE=OFF -DDBG_THREAD_BALANCE=OFF -B build
cd build
make BFS GEN_SRC # BFS_REF
cd ..
# sleep 2

export OMP_NUM_THREADS=$1
export TIMEOUT=60
export NITER=10
export DATASET_DIR="schemas"

# Write source vertices to sources.txt
for dataset in "$DATASET_DIR"/*.json; 
do
  # Remove the .json extension to get the dataset name
  dataset_name=$(basename "$dataset" .json)
  build/GEN_SRC $dataset_name $NITER
done

scripts/run_on_all_datasets.sh $TIMEOUT "runtime_parents_small$OMP_NUM_THREADS" build/BFS small parents
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_parents_large$OMP_NUM_THREADS" build/BFS large parents
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_parents_classic$OMP_NUM_THREADS" build/BFS classic parents
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_distances_small$OMP_NUM_THREADS" build/BFS small distances
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_distances_large$OMP_NUM_THREADS" build/BFS large distances
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_distances_classic$OMP_NUM_THREADS" build/BFS classic distances
scripts/run_on_all_datasets.sh $TIMEOUT "runtime_beamer$OMP_NUM_THREADS" beamer