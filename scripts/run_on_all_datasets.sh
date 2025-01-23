#!/bin/bash

# Check if the executable and number of runs are provided
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "Usage: $0 <num_runs> <timeout> <log_filename> <executable> <executable_args>"
    exit 1
fi

# Function to handle termination signals
terminate_script() {
    echo "Script terminated prematurely."
    exit 1
}

# Trap termination signals
trap terminate_script SIGINT SIGTERM

# TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

NUM_RUNS=$1
EXECUTABLE=$2
TIMEOUT=$3
LOG_FILE=$4 # "output_${TIMESTAMP}.log"
DATASET_DIR="schemas"

# Clear the log file
> $LOG_FILE

# Iterate over all json files in the dataset directory
for dataset in "$DATASET_DIR"/*.json; 
do
    # Remove the .json extension to get the dataset name
    dataset_name=$(basename "$dataset" .json)
    
    # Run the executable with the dataset and append the output to the log file
    EXECUTABLE_ARGS="${@:5}"
    for ((i=1; i<=NUM_RUNS; i++)); do
        echo "Running $EXECUTABLE on $dataset_name with arguments $EXECUTABLE_ARGS (Run $i/$NUM_RUNS)"
        echo "Running $EXECUTABLE on $dataset_name with arguments $EXECUTABLE_ARGS (Run $i/$NUM_RUNS)" >> $LOG_FILE
        echo "Command: $EXECUTABLE $dataset_name $EXECUTABLE_ARGS"
        
        timeout $TIMEOUT $EXECUTABLE "$dataset_name" $EXECUTABLE_ARGS >> $LOG_FILE 2>&1
        if [ $? -eq 124 ]; then
            echo "TIMEOUT"
            echo "TIMEOUT" >> $LOG_FILE
            break
        fi
    done
    echo -e "\n" >> $LOG_FILE
    echo "--------------------------------------------------"
done

echo "All datasets processed. Output logged to $LOG_FILE."