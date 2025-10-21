# Optimizing Breadth-First Search on Modern Energy-Efficient Multicore CPUs

This project contains two implementations of Breadth-First Search (BFS): an OpenMP implementation and a pthreads implementation. Both implementations focus on optimizing BFS performance on modern energy-efficient multicore CPUs by addressing memory access patterns and synchronization overheads.

[Thesis](thesis.pdf) | [Slides](slides.pdf)

(If you are looking for the version presented at the SpeedCode FastCode Challenge 2025 ([paper](https://doi.org/10.1145/3711708.3723452)), please check out the `speedcode-challenge` branch.)
## Overview
Breadth-First Search (BFS) is a fundamental algorithm for graph analysis, but its performance on modern multicore CPUs is generally limited by irregular memory access.

The OpenMP implementation uses an optimized data structure to store graph data, called MergedCSR. By co-locating vertex metadata with its adjacency list, this format improves spatial locality and reduces cache misses.The pthreads implementation uses the MergedCSR structure and some custom components to further improve performance. Namely, it includes a custom implementation of a thread pool, a chunk-based frontier with dynamic work-stealing for load balancing, and a lightweight, custom barrier for scalable synchronization.

Both implementations were evaluated against the GAP Benchmark Suite (included as a submodule) on a diverse set of graphs across x86, RISC-V, and ARM platforms. The results demonstrate that the MergedCSR data structure significantly improves memory performance, enabling a geomean speedup of up to 1.5x over the baseline on large-diameter graphs. Furthermore, the explicit pthreads implementation exhibits superior scalability due to its custom synchronization and outperforms the GAP benchmark with a geomean of 2.28x on road networks and 1.87x on random geometric graphs.

## Requirements
*   A C++17 compatible compiler (e.g., GCC, Clang).
*   [CMake](https://cmake.org/) (version 3.10 or newer).
*   [Make](https://www.gnu.org/software/make/).
*   OpenMP for parallel execution.

The project has been tested with g++ 13.0.0 and OpenMP 4.5.

## Installation

1.  Clone the repository:
    ```sh
    git clone <your-repository-url>
    cd <repository-name>
    ```

2.  Initialize and update the submodules (distributed_mmio and  gapbs):
    ```sh
    git submodule update --init --recursive
    ```

## SbatchMan and MtxMan
This project is compatible with SbatchMan and MtxMan for job scheduling and matrix management. The SbatchMan configuration files are provided in the `scripts` directory. The matrix files used for the evaluation be found in the `scripts/matrices.yaml` directory.

## Running

### Pthreads

The Pthreads implementation can be run from the root directory.

```sh
./pthreads/bin/bfs -f <graph-file.mtx> -n <number-of-runs>
```

*   `-f`: Path to the input graph file in Matrix Market (`.mtx`) format.
*   `-n`: Number of BFS runs to execute.
*   `-s`: Specify a source vertex ID. If not provided, a random source is chosen.

### OpenMP

The OpenMP implementation can also be run from the root directory.

```sh
OMP_NUM_THREADS=<threads> ./openmp/bin/bfs <graph-file.mtx> <runs> <implementation>
```

*   `OMP_NUM_THREADS`: Set the number of OpenMP threads.
*   `<graph-file.mtx>`: Path to the input graph file.
*   `<runs>`: Number of BFS runs.
*   `<implementation>`: The BFS algorithm to use. Options are `reference`, `merged_csr_distances`, `merged_csr_parents`.

### GAP Benchmark Suite (GAPBS)

The GAPBS implementation is located in the gapbs directory.

```sh
cd gapbs
OMP_NUM_THREADS=<threads> ./bfs -f <graph-file> -n <iterations>
```

*   `OMP_NUM_THREADS`: Set the number of OpenMP threads.
*   `-f`: Path to the input graph file. GAPBS supports various formats like `.el`, `.mtx`, `.gr`, and its own serialized `.sg` format.
*   `-n`: Number of iterations to run.