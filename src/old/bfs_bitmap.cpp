#include <graph.hpp>
#include <vector>
#include <profiling.hpp>

namespace benchmark {

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  bool *visited;
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {}
  ~Graph() {
    // destructor logic.
    // If you perform any memory allocations with malloc, new, etc. you must
    // free
    //   them here to avoid memory leaks.
  }

  void BFS(vidType source, weight_type *distances) {
    std::vector<vidType> this_frontier;
    distances[source] = 0;
    visited[source] = true;
    this_frontier.push_back(source);
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto &src : this_frontier) {
        for (uint64_t i = rowptr[src]; i < rowptr[src + 1]; i++) {
          vidType dst = col[i];
          if (!visited[dst]) {
            distances[dst] = distances[src] + 1;
            visited[dst] = true;
            next_frontier.push_back(dst);
          }
        }
      }
      std::swap(this_frontier, next_frontier);
    }
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  Graph* g = new Graph(rowptr, col, N, M);
  g->visited = new bool[N];
  for (uint64_t i = 0; i < N; i++) {
    g->visited[i] = false;
  }
  return g;
}

}