## How to build
To build the project, you need cmake (>=3.10) and a C++17 compiler. Clone the repository and run the following commands in the project's root directory:
```bash
cmake -B build
```
## How to run
To run the project, you need to build it first. After building the project, run the following command in the project's root directory:
```bash
cd build
make <target> SOURCE=<source node>
```
The make target selects the dataset to use. The available targets are listed in the Datasets section. The `ID` of the source node is passed to the program using the `SOURCE` environment variable.

The output is written to a file named `*_distances.out` in the directory where the program is executed from.

## Datasets
| Make target | Dataset           | `\|V\|` | `\|E\|` | Notes          | Source    |
|-------------|-------------------|---------|---------|----------------|-----------|
| test        | test              | 11      | 26      | Handcrafted    |           |
| random      | random_1k_5k      | 1000    | 9940    | Small diameter | Speedcode |
| powergrid   | USpowergrid_n4941 | 4942    | 13190   | Small diameter | [^1]      |

[^1]: https://toreopsahl.com/datasets/#uspowergrid

## Graphs
| ![Test graph](docs/test_graph.png) | 
|:--:| 
| *Test graph* |

