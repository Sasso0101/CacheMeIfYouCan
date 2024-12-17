#include <iostream>
#include "utils.hpp"
#include <parlay/sequence.h>

namespace bfs {

#define MARKED 1 << 31
#define MSB 31

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {}
  ~Graph() {}

  void compute_distances(weight_type *distances, vidType source) {
    for (uint64_t i = 0; i < N; i++) {
      // Check if node has only one neighbour
      if (rowptr[i] == rowptr[i + 1] - 1) {
        vidType neighbour = col[rowptr[i]] & ~(MARKED);
        distances[i] = (col[neighbour] & ~(MARKED)) + 1;
      } else {
        distances[i] = col[rowptr[i]] & ~(MARKED);
      }
    }
    distances[source] = 0;
  }

  void print_frontier(parlay::sequence<vidType> &frontier) {
    std::cout << "Frontier: ";
    for (const auto &v : frontier) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  void BFS(vidType source, weight_type *distances) {
    parlay::sequence<vidType> this_frontier;
    int distance = 0;
    // If the source node has only one neighbor, start the BFS from the neighbor
    if (rowptr[source] == rowptr[source + 1] - 1) {
      this_frontier.push_back(col[rowptr[source]] & ~(MARKED));
      distance = 1;
    } else {
      this_frontier.push_back(rowptr[source]);
    }
    while (!this_frontier.empty()) {
      parlay::sequence<vidType> next_frontier;
      print_frontier(this_frontier);
      for (vidType v = 0; v < this_frontier.size(); ++v) {
        vidType neighbor_index = this_frontier[v];
        // Repeat until all neighbors have been visited
        do {
          vidType neighbor = col[neighbor_index] & ~(MARKED);
          std::cout << "(v: " << this_frontier[v] << ") neighbor: " << neighbor;
          if ((col[neighbor] >> MSB) == 0) {
            col[neighbor] = col[neighbor] | MARKED;
            next_frontier.push_back(neighbor);
            std::cout << " added" << std::endl;
          } else {
            std::cout << " not added" << std::endl;
          }
          neighbor_index++;
        } while ((this_frontier[v] == neighbor_index - 1) ||
                 (col[neighbor_index - 1] >> MSB) == 0);
        col[this_frontier[v]] = distance | MARKED;
      }
      distance++;
      std::swap(this_frontier, next_frontier);
    }
    compute_distances(distances, source);
    std::cout << "Finished BFS!\n";
  }
};

inline void merged_csr(eidType *rowptr, vidType *col, uint64_t N, uint64_t M) {
  for (vidType i = 0; i < N; i++) {
    // Replace list of neighbors with list of neighbor indices in col
    for (uint64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
      col[j] = vidType(rowptr[col[j]]);
    }
    // Mark last neighbor of node
    col[rowptr[i + 1] - 1] = col[rowptr[i + 1] - 1] | MARKED;
  }
}

inline BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                                   uint64_t M) {
  merged_csr(rowptr, col, N, M);
  std::cout << "Merged CSR!\n";
  return new Graph(rowptr, col, N, M);
}

} // namespace bfs