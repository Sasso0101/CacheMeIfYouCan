#include <graph.hpp>
#include <vector>
#include <profiling.hpp>
#include <string>
#ifdef DBG_FRONTIER_SIZE
#include <cstdio>
#endif

namespace reference {

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {}
  ~Graph() {
    // destructor logic.
    // If you perform any memory allocations with malloc, new, etc. you must
    // free
    //   them here to avoid memory leaks.
  }

  void BFS(vidType source, weight_type *distances) {
    #ifdef DBG_FRONTIER_SIZE
      std::vector<vidType> frontier_sizes;
      std::vector<vidType> frontier_max_deg_diff;
      eidType frontier_max_deg;
      eidType frontier_min_deg;
    #endif
    std::vector<vidType> this_frontier;
    distances[source] = 0;
    this_frontier.push_back(source);
    while (!this_frontier.empty()) {
      #ifdef DBG_FRONTIER_SIZE
        frontier_sizes.push_back(this_frontier.size());
        frontier_max_deg = 0;
        frontier_min_deg = std::numeric_limits<eidType>::max();
      #endif
      std::vector<vidType> next_frontier;
      for (const auto &src : this_frontier) {
        #ifdef DBG_FRONTIER_SIZE
          eidType deg = rowptr[src + 1] - rowptr[src];
          if (deg > frontier_max_deg) frontier_max_deg = deg;
          if (deg < frontier_min_deg) frontier_min_deg = deg;
        #endif
        for (uint64_t i = rowptr[src]; i < rowptr[src + 1]; i++) {
          vidType dst = col[i];
          if (distances[src] + 1 < distances[dst]) {
            distances[dst] = distances[src] + 1;
            next_frontier.push_back(dst);
          }
        }
      }
      #ifdef DBG_FRONTIER_SIZE
        frontier_max_deg_diff.push_back(frontier_max_deg - frontier_min_deg);
      #endif
      std::swap(this_frontier, next_frontier);
    }
    #ifdef DBG_FRONTIER_SIZE
      printf("Frontier sizes: ");
      for (const vidType &size : frontier_sizes) {
        printf("%u ", size);
      }
      printf("\n");
      printf("Frontier max deg diff: ");
      for (const eidType &diff : frontier_max_deg_diff) {
        printf("%lu ", diff);
      }
      printf("\n");
    #endif
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M, std::string algorithm) {
  return new Graph(rowptr, col, N, M);
}

}