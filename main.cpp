#include <fstream>
#include <iostream>
#include "bfs.hpp"

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <graph file> <num vertices> <num edges>" << std::endl;
    exit(1);
  }
  // Read graph from file
  std::fstream file;
  file.open("datasets/" + std::string(argv[1]), std::ios::in);
  if (!file) {
    std::cerr << "Unable to open file." << std::endl;
    exit(1);  // terminate with error
  }
  // Load graph from file. The file has format i j, where i and j are the nodes of the edge
  const u_int64_t N = std::stoi(argv[2]);
  const u_int64_t M = std::stoi(argv[3]);
  eidType* rowptr = new eidType[N+1]{};
  vidType* col = new vidType[M]{};
  eidType i;
  vidType j;
  vidType free = 0;
  while (file >> i >> j) {
    col[free] = j;
    free += 1;
    rowptr[i+1] = free;
  }

  // print_graph(rowptr, col, N, M);
  initialize_graph(rowptr, col, N, M);
  print_graph(rowptr, col, N, M);
}