#include "graph.hpp"
#include <benchmark.hpp>
#include <profiling.hpp>
#include <bfs.hpp>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>

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

bool check_correctness(weight_type *distances_ref, weight_type *distances_bfs,
                       uint64_t N) {
  bool correct = true;
  for (uint64_t i = 0; i < N; i++) {
    if (distances_ref[i] != distances_bfs[i]) {
      std::cerr << "Mismatch at node " << i << " ref: " << distances_ref[i]
                << " bfs: " << distances_bfs[i] << std::endl;
      correct = false;
    }
  }
  return correct;
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " <graph file> <num vertices> <num edges> <start node>"
              << std::endl;
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
  vidType start_node = std::stoi(argv[4]);
  vidType free = 0;
  file >> i >> j;
  for (uint64_t k = 0; k < N; k++) {
    while (i == k && file) {
      col[free] = j;
      free += 1;
      file >> i >> j;
    }
    rowptr[k + 1] = free;
  }

  std::cout << "Graph loaded!" << std::endl;

  BaseGraph *g;
  weight_type *distances_ref = new weight_type[N]{};
  for (uint64_t i = 0; i < N; i++) {
    distances_ref[i] = std::numeric_limits<weight_type>::max();
  }
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_THREADINIT;
  auto t1 = std::chrono::high_resolution_clock::now();
  g = benchmark::initialize_graph(rowptr, col, N, M);
  g->BFS(start_node, distances_ref);
  auto t2 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Reference: " << ms_double.count() << "ms\n";

  weight_type *distances_bfs = new weight_type[N]{};
  for (uint64_t i = 0; i < N; i++) {
    distances_bfs[i] = std::numeric_limits<weight_type>::max();
  }
  g = bfs::initialize_graph(rowptr, col, N, M);
  g->BFS(start_node, distances_bfs);

  if (check_correctness(distances_ref, distances_bfs, N)) {
    return 0;
  } else {
    return 1;
  }
  // write_distances(distances_ref, N, "ref_distances.txt");
  // write_distances(distances_bfs, N, "bfs_distances.txt");

  // print_graph(rowptr, col, N, M);
}