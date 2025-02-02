#include <graph.hpp>
#include <iostream>
#include <sstream>

namespace export_graph {

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
    std::ostringstream str;
    for (vidType i = 0; i < N; i++) {
      for (vidType j = rowptr[i]; j < rowptr[i+1]; j++) {
        str << i << " " << col[j] << ", ";
      }
    }
    std::cout << str.str() << std::endl;
  }
};
}