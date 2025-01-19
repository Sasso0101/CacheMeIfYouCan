#!/bin/bash
# Job name
#PBS -N BFS
# Output files
#PBS -o ./stdout.o
#PBS -e ./stderr.e
# Queue name
#PBS -q short_cpuQ
# Set the maximum wall time
#PBS -l walltime=0:06:00
# Number of nodes, cores and memory
#PBS -l select=1:ncpus=96:mem=96gb

cd /home/salvatore.andaloro/Speedcode

# Define the datasets directory
SCHEMA_DIR="./schemas"

# Define the programs to run
PROGRAMS=("merged")

# Define the profiling groups
LIKWID_GROUPS=("L2" "DATA" "L2CACHE" "TLB_DATA" "CYCLE_ACTIVITY" )

OMP_NUM_THREADS=24

# Load modules
module load gcc91 likwid-4.3.4 cmake-3.15.4

rm -rf output
mkdir output

# Loop over each program
for program in "${PROGRAMS[@]}"; do
    # Loop over each dataset in the datasets directory
    for dataset in "$SCHEMA_DIR"/*; do
        # Get the dataset name
        dataset_name=$(basename "$dataset")
        # Remove extension from dataset name
        dataset_name="${dataset_name%.*}"
        
        # Loop over each profiling group
        for likwid_group in "${LIKWID_GROUPS[@]}"; do
            # Run the program with likwid profiling
            echo "${program} with dataset ${dataset_name} and profiling group ${likwid_group}..."
            likwid-perfctr -o "output/${program}-${dataset_name}-${likwid_group}.csv" -C S0 -m -g "$likwid_group" ./build/BFS "$dataset_name" "$program"
        done
    done
done