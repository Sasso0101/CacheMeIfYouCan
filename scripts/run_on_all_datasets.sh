#!/bin/bash

# Check if the executable is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <executable> <executable_args>"
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
    EXECUTABLE_ARGS="${@:2}"
    echo "Running $EXECUTABLE on $dataset_name with arguments $EXECUTABLE_ARGS"
    echo "Running $EXECUTABLE on $dataset_name with arguments $EXECUTABLE_ARGS" >> $LOG_FILE
    echo "Command: $EXECUTABLE $dataset_name $EXECUTABLE_ARGS"
    
    timeout 30 $EXECUTABLE "$dataset_name" $EXECUTABLE_ARGS >> $LOG_FILE 2>&1
    if [ $? -eq 124 ]; then
        echo "TIMEOUT" >> $LOG_FILE
    fi
    echo -e "\n" >> $LOG_FILE
    echo "--------------------------------------------------"
done

echo "All datasets processed. Output logged to $LOG_FILE."