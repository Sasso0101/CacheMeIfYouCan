# Cache-optimized BFS on multi-core CPUs

Breadth-First Search (BFS) performance on shared-memory systems is often limited by irregular memory access and cache inefficiencies. This work presents two optimizations for BFS graph traversal: a bitmap-based algorithm designed for small-diameter graphs and MergedCSR, a graph storage format that improves cache locality for large-scale graphs. Experimental results on real-world datasets show an average 1.3× speedup over a state-of-the-art implementation, with MergedCSR reducing RAM accesses by approximately 15%.

This repository contains the code of the paper "Cache-optimized BFS on multi-core CPUs" by Salvatore D. Andaloro, Thomas Pasquali and Flavio Vella. The paper is available at [TBA]().

## Requirements
The project requires a C++17 compiler. The project has been tested with g++ 12.2.0 and OpenMP 4.5. The project also requires the [CMake](https://cmake.org/) build system. The project has been tested with CMake 3.31.6.

## Building
Clone the repository and run the following commands in the project's root directory:
```bash
cmake -B build
cmake --build build
```

Before running the project, the datasets must be downloaded. This can be done by running the following command in the project's root directory:
```bash
./datasets/download_datasets.sh
```

## Running
To run the project, run the following command in the project's root directory:
```bash
./build/BFS <schema> <source> <implementation> <check>
```
Arguments:
  | Argument   | Description |
  |------------|-----------------------------------------------------------------------------|
  | `<schema>` | Filename of the dataset schema. See the [Datasets](#datasets) section for more details about the available datasets. |
  | `<source>` | Integer. Source vertex of the BFS (`0` by default) |
  | `<algorithm>` | Implementation used to perform the BFS. One of `merged_csr_parents`, `merged_csr`, `bitmap`, `classic`, `reference` or `heuristic` (`heuristic` by default). See the paper for more details on the implementations. |
  | `<check>`  | `true` or `false`. Checks correctness of the result using a simple single-threaded implementation. (`false` by default) |

## Testing

To run the tests, run the following command in the project's root directory:
```bash
ctest --test-dir build/tests
```
Each test is a BFS run on a small dataset with a different implementation. The test checks the correctness of the result by comparing it with the reference implementation.

## Datasets
The datasets used are provided by the [PPoPP'25 FastCode Challenge](https://fastcode.org/events/fastcode-challenge/spe4ic/#dataset-diversity). Note that they weight several gigabytes in total, so they may take a while to download.
|        Graph Name       | \|V\| |  \|E\| |      Notes     | Filename |
|-------------------------|:-----:|:------:|:--------------:|----------|
| Social_Network_1        | 4.9M  | 85.8M  | Small diameter | Social_Network_1.json |
| Web_Graph_1             | 6.6M  | 300M   | Small diameter | Web_Graph_1.json |
| Collaboration_Network_1 | 1.1M  | 113M   | Small diameter | Collaboration_Network_1.json |
| Synthetic_Dense_1       | 10M   | 1B     | Small diameter | Synthetic_Dense_1.json |
| Road_Network_1          | 22.1M | 30.0M  | Large diameter | Road_Network_1.json |
| Road_Network_2          | 87.0M | 112.9M | Large diameter | Road_Network_2.json |
| kNN_Graph_1             | 24.9M | 158M   | Large diameter | kNN_Graph_1.json |
| Synthetic_Sparse_1      | 10M   | 40M    | Large diameter | Synthetic_Sparse_1.json |

#### Acknowledgements
This work was partially supported by the EuroHPC JU project within Net4Exa project under grant agreement No 101175702.