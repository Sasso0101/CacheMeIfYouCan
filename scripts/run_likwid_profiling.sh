#!/bin/bash
# Job name
#PBS -N BFS_LIKWID_PROFILING
# Output files
#PBS -o ./likwid_stdout.o
#PBS -e ./likwid_stderr.e
# Queue name
#PBS -q short_cpuQ
# Set the maximum wall time
#PBS -l walltime=06:00:00
# Number of nodes, cores and memory
#PBS -l select=1:ncpus=96:mem=100gb
#PBS -l place=exclhost

# Function to handle termination signals
terminate_script() {
    echo "Script terminated prematurely."
    exit 1
}

# Trap termination signals
trap terminate_script SIGINT SIGTERM

cd /home/salvatore.andaloro/CacheMeIfYouCan

# Define the datasets directory
SCHEMA_DIR="./schemas"

# Define the programs to run
PROGRAMS=("classic" "large" "small")

# Define the profiling groups
# LIKWID_GROUPS=("L2" "DATA" "L2CACHE" "TLB_DATA" "CYCLE_ACTIVITY")
LIKWID_GROUPS=("L2CACHE" "L3CACHE")

OMP_NUM_THREADS=24

# Load modules
module load gcc91 likwid-4.3.4 cmake-3.15.4

rm -rf output
mkdir output

# Loop over each program
for program in "${PROGRAMS[@]}"; do
    for threads in 24 48 96; do
        export OMP_NUM_THREADS=$threads
        # Loop over each dataset in the datasets directory
        for dataset in "$SCHEMA_DIR"/*; do
            # Get the dataset name
            dataset_name=$(basename "$dataset")
            # Remove extension from dataset name
            dataset_name="${dataset_name%.*}"
            
            # Loop over each profiling group
            for likwid_group in "${LIKWID_GROUPS[@]}"; do
                # Run the program with likwid profiling
                echo "Running BFS ${program} (${threads} threads) on dataset ${dataset_name} and profiling group ${likwid_group}..."
                echo "likwid-perfctr -o output/${program}-${threads}-${dataset_name}-${likwid_group}.csv -C 0-$((OMP_NUM_THREADS-1)) -m -g $likwid_group ./build/BFS $dataset_name $program"
                likwid-perfctr -o "output/${program}-${threads}-${dataset_name}-${likwid_group}.csv" -C "0-$((OMP_NUM_THREADS-1))" -m -g "$likwid_group" ./build/BFS "$dataset_name" "$program"
            done
        done
    done
done