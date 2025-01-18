#include <cassert>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <profiling.hpp>

namespace test_openmp {

class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {
  }
  ~Graph() {}

  void BFS(vidType source, weight_type *distances) {
    // iteratively remove leaves from graph
    int leaves = 0;
    int *potential_leaves = new int[N];
    // init potential leaves to 0
    for (uint64_t i = 0; i < N; i++) {
      potential_leaves[i] = 0;
    }
    for (uint64_t i = 0; i < N; i++) {
      if (rowptr[i + 1] - rowptr[i] == 1) {
        vidType neighbor = col[rowptr[i]];
        // find leaf in parent and mark it
        if (neighbor != -1) {
          vidType to_clear = rowptr[neighbor];
          while (to_clear < rowptr[neighbor + 1]) {
            if (to_clear != -1 && col[to_clear] == i) break;
            to_clear++;
          }
          if (to_clear < rowptr[neighbor + 1]) {
            col[to_clear] = -1;
            leaves++;
            potential_leaves[neighbor] += 1;
          }
        } else {
          potential_leaves[i] = 0;
          std::cout << "2 connected leaves " << i << std::endl;
        }
      }
    }
    std::cout << "Level 1 leaves: " << leaves << std::endl;
    for (int level = 2; level < 4; level++) {
      leaves = 0;
      // compute potential leaves
      for (uint64_t i = 0; i < N; i++) {
        if ((rowptr[i+1] - rowptr[i]) > 1 && potential_leaves[i] == (rowptr[i+1] - rowptr[i] - 1)) {
          leaves++;
          // find non-leaf in list
          vidType parent = rowptr[i];
          while (col[parent] == -1) {
            parent++;
          }
          vidType neighbor = col[parent];
          vidType to_clear = rowptr[neighbor];
          while (col[to_clear] != i) {
            to_clear++;
          }
          col[to_clear] = -1;
          potential_leaves[neighbor] += 1;
          potential_leaves[i] = 0;
        } else if ((rowptr[i+1] - rowptr[i]) > 1 && potential_leaves[i] == (rowptr[i+1] - rowptr[i])) {
          leaves++;
          potential_leaves[i] = 0;
          std::cout << "root " << i << std::endl;
        }
      }
      std::cout << "Level " << level << " leaves: " << leaves << std::endl;
    }
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  return new Graph(rowptr, col, N, M);
}

}