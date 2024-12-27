#include <chrono>
#include <iostream>
#include <graph.hpp>
#include <parlay/sequence.h>
#include <omp.h>

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

  inline vidType copy_unmarked(vidType i) {
    return col[i] & ~(MARKED);
  }

  inline bool is_marked(vidType i) {
    return (col[i] >> MSB) != 0;
  }

  inline void mark(vidType i) {
    col[i] = col[i] | MARKED;
  }

  inline void set_distance(vidType i, weight_type distance) {
    col[i] = distance | MARKED;
  }

  inline bool is_leaf(vidType i) {
    return (rowptr[i] == rowptr[i + 1] - 1);
  }

  void compute_distances(weight_type *distances, vidType source) {
    for (uint64_t i = 0; i < N; i++) {
      // Check if node has only one neighbor
      if (is_leaf(i)) {
        vidType neighbour = copy_unmarked(rowptr[i]);
        if (is_marked(neighbour)) {
          distances[i] = copy_unmarked(neighbour) + 1;
        }
      } else if (is_marked(rowptr[i])) {
        distances[i] = copy_unmarked(rowptr[i]);
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
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<vidType> this_frontier;
    weight_type distance = 0;
    // If the source node has only one neighbor, start the BFS from the neighbor
    if (rowptr[source] == rowptr[source + 1] - 1) {
      this_frontier.push_back(copy_unmarked(rowptr[source]));
      distance = 1;
    } else {
      this_frontier.push_back(rowptr[source]);
    }
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (vidType v = 0; v < this_frontier.size(); ++v) {
        vidType neighbor_index = this_frontier[v];
        // Repeat until all neighbors have been visited
        do {
          vidType neighbor = copy_unmarked(neighbor_index);
          if (!is_marked(neighbor)) {
            mark(neighbor);
            next_frontier.push_back(neighbor);
          }
          neighbor_index++;
        } while ((this_frontier[v] == neighbor_index - 1) ||
                 !is_marked(neighbor_index - 1));
        set_distance(this_frontier[v], distance);
      }
      distance++;
      std::swap(this_frontier, next_frontier);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "BFS: " << ms_double.count() << "ms\n";

    t1 = std::chrono::high_resolution_clock::now();
    compute_distances(distances, source);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Postprocessing: " << ms_double.count() << "ms\n";
  }
};

void merged_csr(eidType *rowptr, vidType *col, uint64_t N, uint64_t M) {
  for (vidType i = 0; i < N; i++) {
    // Replace list of neighbors with list of neighbor indices in col
    for (uint64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
      col[j] = vidType(rowptr[col[j]]);
    }
    // Mark last neighbor of node
    col[rowptr[i + 1] - 1] = col[rowptr[i + 1] - 1] | MARKED;
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                                   uint64_t M) {
  auto t1 = std::chrono::high_resolution_clock::now();
  merged_csr(rowptr, col, N, M);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Preprocessing: " << ms_double.count() << "ms\n";
  return new Graph(rowptr, col, N, M);
}

}