#!/bin/bash

# Check if the executable is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <executable>"
    exit 1
fi

EXECUTABLE=$1
DATASET_DIR="schemas"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="output_$EXECUTABLE_$TIMESTAMP.log"

# Clear the log file
> $LOG_FILE

# Iterate over all json files in the dataset directory
for dataset in "$DATASET_DIR"/*.json; 
do
    # Remove the .json extension to get the dataset name
    dataset_name=$(basename "$dataset" .json)
    
    # Run the executable with the dataset and append the output to the log file
    echo "Running $EXECUTABLE on $dataset_name"
    echo "Running $EXECUTABLE on $dataset_name" >> $LOG_FILE
    echo "Command: $EXECUTABLE $dataset_name"
    $EXECUTABLE "$dataset_name" >> $LOG_FILE 2>&1
    echo -e "\n" >> $LOG_FILE
    echo "--------------------------------------------------"
done

echo "All datasets processed. Output logged to $LOG_FILE."