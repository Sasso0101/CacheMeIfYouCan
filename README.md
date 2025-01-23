## How to build
To build the project, you need cmake (>=3.10) and gcc9.1.0. Clone the repository and run the following commands in the project's root directory:
```bash
cmake -B build
cd build
make BFS
```

## How to run
To run the project, run the following command in the project's root directory:
```bash
./build/BFS <dataset> <algorithm>
```
The available datasets are listed in the [Datasets](#datasets) section.

The available algorithms are `merged` and `nomerged`. The `merged` algorithm merges the list of vertices and the list of edges into a array, while the `nomerged` algorithm keeps them separate, but uses a bitmap to mark visited vertices.

## Datasets
The datasets are automatically downloaded in the `datasets` directory when the project is built. Note that they weigh several gigabytes in total, so they may take a while to download. The datasets are listed in the table below.
|        Graph Name       | \|V\| |  \|E\| |      Notes     |
|:-----------------------:|:-----:|:------:|:--------------:|
| Social_Network_1        | 4.9M  | 85.8M  | Small diameter |
| Web_Graph_1             | 6.6M  | 300M   | Small diameter |
| Collaboration_Network_1 | 1.1M  | 113M   | Small diameter |
| Synthetic_Dense_1       | 10M   | 1B     | Small diameter |
| Road_Network_1          | 22.1M | 30.0M  | Large diameter |
| Road_Network_2          | 87.0M | 112.9M | Large diameter |
| kNN_Graph_1             | 24.9M | 158M   | Large diameter |
| Synthetic_Sparse_1      | 10M   | 40M    | Large diameter |

## Profiling
The project is already set up to use the [LIKWID](https://github.com/RRZE-HPC/likwid) suite for profiling. The LIKWID suite must be installed on the system. To install the LIKWID suite, follow the instructions in the [README file](https://github.com/RRZE-HPC/likwid?tab=readme-ov-file#download-build-and-install) in the LIKWID repository. To list the available profiling groups, run `likwid-perfctr -a`. To view the detailed description of a group, run `likwid-perfctr -g <group> -H`.

For example, to view the `CYCLE_ACTIVITY` statistics (measures cycles spent waiting for data from the cache and memory hierarchy), run the following command:
```bash
cd build
likwid-perfctr -C 0 -g CYCLE_ACTIVITY -m ./BFS Collaboration_Network_1 merged
```


## Python env setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib itables pandas
```

## CMake flags
```bash
cmake -DDBG_FRONTIER_SIZE=OFF -DDBG_THREAD_BALANCE=ON ..
```