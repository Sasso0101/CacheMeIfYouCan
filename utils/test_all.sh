#!/bin/bash
cd ..
for i in {0..999}
do
  build/BFS datasets/random_1k_5k.txt 1000 9940 $i
  if [ $? -ne 0 ]; then
    echo "Test failed for $i"
    exit 1
  fi
done