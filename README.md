# Cache-optimized BFS on multi-core CPUs

Breadth-First Search (BFS) performance on shared-memory systems is often limited by irregular memory access and cache inefficiencies. This work presents two optimizations for BFS graph traversal: a bitmap-based algorithm designed for small-diameter graphs and MergedCSR, a graph storage format that improves cache locality for large-scale graphs. Experimental results on real-world datasets show an average 1.3× speedup over a state-of-the-art implementation, with MergedCSR reducing RAM accesses by approximately 15%.

This repository contains the code of the paper "Cache-optimized BFS on multi-core CPUs" by Salvatore D. Andaloro, Thomas Pasquali and Flavio Vella. The paper is available at [TBA]().

## Requirements
The project requires a C++17 compiler. The project has been tested with g++ 12.2.0 and OpenMP 4.5. The project also requires the [CMake](https://cmake.org/) build system. The project has been tested with CMake 3.31.6.

## Building
Clone the repository and run the following commands in the project's root directory:
```bash
cmake -B build
make -C build BFS
```
Download the suggested datasets using MtxMan:
```bash
mtxman sync yaml_files/matrices.yaml
```

## Running
To run the project, run the following command in the project's root directory:
```bash
./build/BFS <schema> <source> <implementation> <check>
```
Arguments:
  | Argument   | Description |
  |------------|-----------------------------------------------------------------------------|
  | `<matrix>` | Filename of matrix (in .mtx or .bmtx format). |
  | `<runs>` | Integer. Number of runs to execute (default: 1) |
  | `<algorithm>` | Implementation used to perform the BFS. One of `merged_csr_parents`, `merged_csr`, `bitmap`, `classic`, `reference` or `heuristic` (`heuristic` by default). See the paper for more details on the implementations. |
  | `<check>`  | `true` or `false`. Checks correctness of the result using a simple single-threaded implementation. (`false` by default) |

## Acknowledgements
This work was partially supported by the EuroHPC JU project within Net4Exa project under grant agreement No 101175702.