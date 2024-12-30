#include <chrono>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <vector>
#include <profiling.hpp>

namespace bfs_cached {

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

  inline vidType copy_unmarked(vidType i) { return rowptr[i] & ~(MARKED); }

  inline bool is_marked(vidType i) { return (rowptr[i] >> MSB) != 0; }

  inline void mark(vidType i) { rowptr[i] = rowptr[i] | MARKED; }

  inline void set_distance(vidType i, weight_type distance) {
    rowptr[i] = distance | MARKED;
  }

  inline bool is_unconnected(vidType i) { return (rowptr[i] == rowptr[i + 1]); }

  inline bool is_leaf(vidType i) { return (rowptr[i] == rowptr[i + 1] - 1); }

  void compute_distances(weight_type *distances, vidType source) {
    for (uint64_t i = 0; i < N; i++) {
      if (is_marked(i)) {
        distances[i] = copy_unmarked(i);
      }
    }
    distances[source] = 0;
  }

  void print_frontier(std::vector<vidType> &frontier) {
    std::cout << "Frontier: ";
    for (const auto &v : frontier) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  void BFS(vidType source, weight_type *distances) {
    auto t1 = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("BFS");
    std::vector<vidType> this_frontier;
    weight_type distance = 0;
    this_frontier.push_back(source);

    if (is_unconnected(source)) {
      distances[source] = 0;
      return;
    }
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto &src : this_frontier) {
        uint64_t i = rowptr[src] & ~(MARKED);
        vidType dst;
        for (; (col[i] >> MSB) == 0; i++) {
          dst = col[i];
          if (!is_marked(dst)) {
            mark(dst);
            next_frontier.push_back(dst);
          }
        }
        dst = col[i] & ~(MARKED);
        if (!is_marked(dst)) {
          mark(dst);
          next_frontier.push_back(dst);
        }
        set_distance(src, distance);
      }
      distance++;
      std::swap(this_frontier, next_frontier);
    }
    LIKWID_MARKER_STOP("BFS");
    LIKWID_MARKER_CLOSE;
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
  vidType start = 0;
  // Skip first rows that could contain zero
  while (rowptr[start+1] == 0) {
    start++;
  }
  for (vidType i = start; i < N; i++) {
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