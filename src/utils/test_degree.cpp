#include <cassert>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <profiling.hpp>

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
    int average = 0;
    int max = 0;
    int min = std::numeric_limits<int>::max();
    int *histogram = new int[100]{};
    // compute histogram of degrees
    for (uint64_t i = 0; i < N; i++) {
      average += rowptr[i + 1] - rowptr[i];
      if (rowptr[i + 1] - rowptr[i] > max) {
        max = rowptr[i + 1] - rowptr[i];
      }
      if (rowptr[i + 1] - rowptr[i] < min) {
        min = rowptr[i + 1] - rowptr[i];
      }
    }
    for (uint64_t i = 0; i < N; i++) {
      histogram[(rowptr[i + 1] - rowptr[i]) / (max / 100)] += 1;
    }
    float average_f = (float)average / N;
    std::cout << "Average degree: " << average_f << std::endl;
    std::cout << "Max degree: " << max << std::endl;
    std::cout << "Min degree: " << min << std::endl;
    for (int i = 0; i < 100; i++) {
      std::cout << i << "," << histogram[i] << std::endl;
    }
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  return new Graph(rowptr, col, N, M);
}