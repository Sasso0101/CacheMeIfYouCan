## How to build
To build the project, you need cmake (>=3.10) and a C++17 compiler. Clone the repository and run the following commands in the project's root directory:
```bash
cd build
cmake .. -DPARLAY_OPENMP=On -DCMAKE_EXPORT_COMPILE_COMMANDS=1
```
## How to run
To run the project, run the following commands in the project's root directory:
```bash
cd build
make <target> SOURCE=<source node>
```
The make target selects the dataset to use. The available targets are listed in the Datasets section. The `ID` of the source node is passed to the program using the `SOURCE` environment variable.

The output is written to a file named `*_distances.out` in the directory where the program is executed from.

## Datasets
| Make target | Filename          | `\|V\|` | `\|E\|`  | Notes            | Source    |
|-------------|-------------------|---------|----------|------------------|-----------|
| test        | test.txt          | 16      | 32       | Handcrafted      |           |
| random      | random_1k_5k.txt  | 1000    | 9940     | Small diameter   | Speedcode |
| powergrid   | powergrid.txt     | 4942    | 13190    | Small diameter   | [^1]      |
| epinions    | epinions.txt      | 75879   | 508837   | Made undirected  | [^2]      |
| twitch      | twitch.txt        | 168113  | 13595114 | Made undirected  | [^3]      |

[^1]: https://toreopsahl.com/datasets/#uspowergrid
[^2]: https://www.kaggle.com/datasets/wolfram77/graphs-social?resource=download
[^3]: https://snap.stanford.edu/data/twitch_gamers.html

## Graphs
| ![Test graph](docs/test_graph.png) | 
|:--:| 
| *Test graph* |

