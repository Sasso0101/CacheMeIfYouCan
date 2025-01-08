## How to build
To build the project, you need cmake (>=3.10) and clang++. Clone the repository and run the following commands in the project's root directory:
```bash
cd build
cmake .. -DPARLAY_OPENMP=On
```
The `-DPARLAY_OPENMP=On` flag enables OpenMP support in the Parlay library.

## How to run
To run the project, run the following commands in the project's root directory:
```bash
cd build
make <target> SOURCE=<source node>
```

The make target selects the dataset to use. The available targets are listed in the Datasets section. The `ID` of the source node is passed to the program using the `SOURCE` environment variable.

Each time the code is compiled, a vectorization report is written to a file named `optinfo.txt` in the build directory.

The program will output the computed distances to a file named `*_distances.out` in the directory where the `make` command is executed from.

## Datasets
The datasets are automatically downloaded in the `datasets` directory when the project is built. Note that the pokec dataset is 600MB in size and may take a while to download.
| Make target | Filename          | `\|V\|` | `\|E\|`  | Notes            | Source    |
|-------------|-------------------|---------|----------|------------------|-----------|
| test        | test.txt          | 16      | 32       |                  | Handcrafted |
| unconnected | unconnected.txt   | 4       | 2        | Contains unconnected vertices | Handcrafted |
| random      | random_1k_5k.txt  | 1000    | 9940     |                  | Speedcode |
| powergrid   | powergrid.txt     | 4942    | 13190    | Power grid       | [^1]      |
| epinions    | epinions.txt      | 75879   | 508837   | Social network   | [^2]      |
| twitch      | twitch.txt        | 168113  | 13595114 | Social network   | [^3]      |
| pokec       | pokec.txt         | 1632804 | 44603928 | Social network   | [^4]      |
| roadnet-ca  | roadnet-ca.txt    | 1971281 | 5533214  | Road network (California) | [^5] |

The datasets taken from the SNAP dataset collection were converted to undirected graphs by adding the reverse of each edge. Then, the edges were sorted by source vertex and destination vertex. The Python script used to convert the datasets is available in the `utils` directory.

[^1]: https://toreopsahl.com/datasets/#uspowergrid
[^2]: http://snap.stanford.edu/data/soc-Epinions1.html
[^3]: https://snap.stanford.edu/data/twitch_gamers.html
[^4]: http://snap.stanford.edu/data/soc-Pokec.html
[^5]: http://snap.stanford.edu/data/roadNet-CA.html

## Graphs
<table width="100%">
  <tbody>
    <tr>
      <td width="50%"><img src="docs/test_graph.png"/></td>
      <td width="50%"><img src="docs/unconnected_graph.png"/></td>
    </tr>
    <tr align="center">
      <td width="50%"><i>Test graph</i></td>
      <td width="50%"><i>Graph with unconnected vertices</i></td>
    </tr>
  </tbody>
</table>

## Profiling
The project is already set up to use the [LIKWID](https://github.com/RRZE-HPC/likwid) suite for profiling. The LIKWID suite must be installed on the system. To install the LIKWID suite, follow the instructions in the [README file](https://github.com/RRZE-HPC/likwid?tab=readme-ov-file#download-build-and-install) in the LIKWID repository. To list the available profiling groups, run `likwid-perfctr -a`. To view the detailed description of a group, run `likwid-perfctr -g <group> -H`.

For example, to view the `CYCLE_ACTIVITY` statistics (measures cycles spent waiting for data from the cache and memory hierarchy), run the following command:
```bash
cd build
likwid-perfctr -C 0 -g CYCLE_ACTIVITY -m ./BFS ../datasets/pokec.txt 1632804 44603928 10
```
