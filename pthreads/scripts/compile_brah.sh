#!/bin/bash

set -e

impls=(atomic main openmp)

for impl in "${impls[@]}"; do

    cd $impl
    implname=$impl
    [[ $implname == "main" ]] && implname="mutex"
    mkdir -p targets

    if [[ $impl == "openmp" ]]; then
        cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -B build
        make -j8 -C build BFS
        mv build/BFS "targets/bfs_${implname}"
        echo "Built: targets/bfs_${implname}"
        continue
    fi

    for chunksize in 4 8 16 32 64 256 1024 4096; do
        for ncpus in 1 2 4 8 16 32 64; do
            # [[ -f "targets/bfs_${implname}_chunksize${chunksize}_${ncpus}cpus" ]] && continue
            make clean
            CHUNK_SIZE=${chunksize} MAX_THREADS=${ncpus} CC=gcc CXX=g++ make -j8 bin/bfs
            mv bin/bfs "targets/bfs_${implname}_chunksize${chunksize}_${ncpus}cpus"
            echo "Built: targets/bfs_${implname}_chunksize${chunksize}_${ncpus}cpus"
        done
    done
    cd ..
done
