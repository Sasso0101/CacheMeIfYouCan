#!/bin/bash

# Check if the executable and number of runs are provided
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <timeout> <log_filename> <executable> <executable_args>"
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

TIMEOUT=$1
LOG_FILE="$2.log" # "output_${TIMESTAMP}.log"
EXECUTABLE=$3

# Read source vertices from sources.txt, the file has the following format:
# <dataset_name> <source_node_1> <source_node_2> ... <source_node_n>
declare -A DATASETS
while IFS= read -r line
do
    # Split the line into an array
    IFS=' ' read -r -a array <<< "$line"
    dataset_name=${array[0]}
    DATASETS[$dataset_name]=${array[@]:1}
done < sources.txt

# Clear the log file
> $LOG_FILE

# Iterate over datasets in sources.txt
for dataset in "${!DATASETS[@]}"; 
do
    # Run the executable with the dataset and append the output to the log file
    EXECUTABLE_ARGS="${@:4}"
    # Enumerate the sources for the dataset
    run_count=1
    length=$(echo ${DATASETS[$dataset]} | wc -w)
    for source in ${DATASETS[$dataset]}; do
        echo "Running $EXECUTABLE on $dataset with arguments $EXECUTABLE_ARGS with source $source (Run $run_count/$length)"
        echo "Running $EXECUTABLE on $dataset with arguments $EXECUTABLE_ARGS with source $source (Run $run_count/$length)" >> $LOG_FILE
        
        if [ "$EXECUTABLE" == "beamer" ]; then
            echo "Command: gapbs/bfs -f datasets/$dataset.sg -n 1 -r $source"
            gapbs/bfs -f datasets/"${dataset}".sg -n 1 -r "$source" >> $LOG_FILE 2>&1
        else
            echo "Command: $EXECUTABLE $dataset $source $EXECUTABLE_ARGS"
            "$EXECUTABLE" "$dataset" "$source" ${EXECUTABLE_ARGS[@]} >> $LOG_FILE 2>&1
        fi
        run_count=$((run_count + 1))
    done
    echo -e "\n" >> $LOG_FILE
    echo "--------------------------------------------------"
done

echo "All datasets processed. Output logged to $LOG_FILE."