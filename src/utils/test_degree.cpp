#include <cassert>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <profiling.hpp>
#include <algorithm>
#include "../input.hpp"

#define HIST_SIZE 100

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
    eidType tot_degree = 0;
    eidType max_deg = 0;
    eidType min_deg = std::numeric_limits<eidType>::max();
    eidType degrees[N];
    eidType bins[HIST_SIZE] = {0};
    eidType frequencies[HIST_SIZE] = {0};
    
    for (uint64_t i = 0; i < N; i++) {
      eidType deg = rowptr[i + 1] - rowptr[i];
      tot_degree += deg;
      max_deg = std::max(max_deg, deg);
      min_deg = std::min(min_deg, deg);
      degrees[i] = deg;
    }

    eidType bin_width = (max_deg - min_deg + HIST_SIZE - 1) / HIST_SIZE; // Ceiling division
    for (int i = 0; i < HIST_SIZE; i++) {
      bins[i] = min_deg + i * bin_width;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < HIST_SIZE; j++) {
            if (degrees[i] >= bins[j] && degrees[i] < bins[j] + bin_width) {
                frequencies[j]++;
                break;
            }
        }
    }
    
    std::cout << "Average degree: " << ((double)tot_degree / N) << std::endl;
    std::cout << "Max degree: " << max_deg << std::endl;
    std::cout << "Min degree: " << min_deg << std::endl;
    std::cout << "Histogram Data: ";
    for (int i = 0; i < HIST_SIZE; i++) {
      printf("(%d-%d:%d)", bins[i], bins[i] + bin_width - 1, frequencies[i]);
    }
    // for (int i = 0; i < N; i++) {
    //   std::cout << degrees[i] << " ";
    // }
    std::cout << std::endl;
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M, std::string s) {
  return new Graph(rowptr, col, N, M);
}

int main(const int argc, char **argv) {
  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, "", "");
  p.run();
}