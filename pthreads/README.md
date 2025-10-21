# Thesis

Optimized BFS implementation in C using pthreads.

## Prerequisites
A Linux system with the following software installed:
- git
- gcc with support for C11
- cmake
- make

## Installation

Install the submodules required for the project:
```bash
git submodule update --init --recursive
```

Run make in the root directory to build the project:
  ```bash
  make
  ```

## Usage

To run the BFS implementation, use the following command:
```bash
./bin/bfs -f <input file> -n <number of runs>
```
For more options, use the `-h` flag:
```bash
./bin/bfs -h
```

Note: the input file must be in the `.mtx` format. You can find many graphs in this format in the [SuiteSparse Matrix Collection](https://sparse.tamu.edu/).