#include <cmath>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>

typedef uint64_t vidType;
typedef uint64_t eidType;
typedef uint64_t weight_type;

class BaseGraph {
 public:
  virtual ~BaseGraph() {}
  virtual void BFS(vidType source, weight_type *distances) = 0;
};

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
    std::vector<vidType> this_frontier;
    distances[source] = 0;
    this_frontier.push_back(source);
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto &src : this_frontier) {
        for (uint64_t i = rowptr[src]; i < rowptr[src + 1]; i++) {
          vidType dst = col[i];
          if (distances[src] + 1 < distances[dst]) {
            distances[dst] = distances[src] + 1;
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
  std::unordered_multiset<eidType> to_check;
  uint64_t removed = 0;
  for (uint64_t i = 0; i < N; i++) {
    // Tree detected
    if (rowptr[i + 1] - rowptr[i] == 1) {
      to_check.insert(col[rowptr[i]]);
      rowptr[i] = -1;
      removed++;
    }
  }

  for (auto v = to_check.begin(); v != to_check.end();) {
    uint64_t count = to_check.count(*v);
    if (rowptr[*v + 1] - rowptr[*v] == count) {
      for (uint64_t i = rowptr[*v]; i < rowptr[*v + 1]; i++) {
        col[rowptr[i]] = -1;
      }
      removed++;
    }
    std::advance(v, count);
  }
  std::cout << "Removed " << removed << " nodes" << std::endl;

  return new Graph(rowptr, col, N, M);
}
