#!/bin/bash

# Define the datasets directory
SCHEMA_DIR="./schema"

# Define the programs to run
PROGRAMS=("merged nomerged")

# Define the profiling groups
GROUPS=("L2" "DATA")

# Loop over each program
for program in "${PROGRAMS[@]}"; do
    # Loop over each dataset in the datasets directory
    for dataset in "$SCHEMA_DIR"/*; do
        # Get the dataset name
        dataset_name=$(basename "$dataset")
        
        # Loop over each profiling group
        for group in "${GROUPS[@]}"; do
            # Run the program with likwid profiling
            likwid-perfctr -o "$program-$dataset_name.csv" -C 0 -g "$group" ./"$program" "$dataset_name"
        done
    done
done