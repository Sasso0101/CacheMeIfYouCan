#!/bin/bash

set -e

module load llvm/cross/EPI-development cmake/3.23.2

impls=(atomic main openmp)

for impl in "${impls[@]}"; do

    cd $impl
    implname=$impl
    [[ $implname == "main" ]] && implname="mutex"
    mkdir -p targets

    if [[ $impl == "openmp" ]]; then
        cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -B build
        make -j16 -C build BFS
        mv build/BFS "targets/bfs_${implname}"
        echo "Built: targets/bfs_${implname}"
        continue
    fi

    for chunksize in 4 8 16 32 64 256 1024 4096; do
        for ncpus in 1 2 4 8 16 32 64; do
            make clean
            CHUNK_SIZE=${chunksize} MAX_THREADS=${ncpus} CC=clang CXX=clang++ make -j16 bin/bfs
            mv bin/bfs "targets/bfs_${implname}_chunksize${chunksize}_${ncpus}cpus"
            echo "Built: targets/bfs_${implname}_chunksize${chunksize}_${ncpus}cpus"
        done
    done
    cd ..
done
