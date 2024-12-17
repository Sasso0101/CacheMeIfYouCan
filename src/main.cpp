#include "benchmark.hpp"
#include "bfs.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>

#define OUTPUT_FOLDER "output/"

void write_col(vidType *col, uint64_t M, std::string filename) {
  std::ofstream file(filename);
  for (uint64_t i = 0; i < M; i++) {
    file << "i: " << i << " " << (col[i] & ~(1 << 31));
    if (col[i] >> 31) {
      file << " x";
    }
    file << std::endl;
  }
}

void print_graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M) {
  for (uint64_t i = 0; i < N; i++) {
    for (uint64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
      std::cout << i << " (j: " << j << ") " << (col[j] & (0 << 31))
                << std::endl;
    }
  }
}

void write_distances(weight_type *distances, uint64_t N, std::string filename) {
  std::ofstream file(filename);
  for (uint64_t i = 0; i < N; i++) {
    file << i << ": " << distances[i] << std::endl;
  }
}

void test_bfs_implementation(BaseGraph *(*initialize_graph)(eidType *,
                                                            vidType *, uint64_t,
                                                            uint64_t),
                             eidType *rowptr, vidType *col, uint64_t N,
                             uint64_t M, const std::string &output_filename, vidType start_node) {
  BaseGraph *g = initialize_graph(rowptr, col, N, M);

  weight_type *distances = new weight_type[N]{};
  for (uint64_t i = 0; i < N; i++) {
    distances[i] = std::numeric_limits<weight_type>::max();
  }
  std::cout << "Start node " << start_node << std::endl;
  write_col(col, M, output_filename+"_col.out");
  g->BFS(start_node, distances);
  write_distances(distances, N, output_filename+"_distances.out");
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " <graph file> <num vertices> <num edges> <start node>" << std::endl;
    exit(1);
  }
  // Read graph from file
  std::fstream file;
  file.open(std::string(argv[1]), std::ios::in);
  if (!file) {
    std::cerr << "Unable to open file." << std::endl;
    exit(1);
  }

  // Load graph from file. The file format is assumed to be:
  // <source vertex> <destination vertex>
  const u_int64_t N = std::stoi(argv[2]);
  const u_int64_t M = std::stoi(argv[3]);
  eidType *rowptr = new eidType[N + 1]{};
  vidType *col = new vidType[M]{};
  eidType i;
  vidType j;
  vidType free = 0;
  while (file >> i >> j) {
    col[free] = j;
    free += 1;
    rowptr[i + 1] = free;
  }

  std::cout << "Graph loaded!" << std::endl;

  test_bfs_implementation(benchmark::initialize_graph, rowptr, col, N, M, "benchmark", std::stoi(argv[4]));
  test_bfs_implementation(bfs::initialize_graph, rowptr, col, N, M, "bfs", std::stoi(argv[4]));

  // print_graph(rowptr, col, N, M);
}