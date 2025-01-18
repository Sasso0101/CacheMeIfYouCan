#include "input.hpp"
#include <complete.hpp>
#include <profiling.hpp>
#include <complete_nomerged.hpp>

int main(const int argc, char **argv) {
  std::ifstream in("schema/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  std::string code_version = std::string(argv[2]);
  BaseGraph *(*init)(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  if (code_version == "nomerged") {
    init = nomerged::initialize_graph;
  } else {
    init = complete::initialize_graph;
  }
  ProblemInput p = ProblemInput(data, init);
  LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
    LIKWID_MARKER_THREADINIT;
  }
  #pragma omp parallel
  {
    LIKWID_MARKER_START("bfs");
  }
  p.run();
  #pragma omp parallel
  {
    LIKWID_MARKER_STOP("bfs");
  }
  LIKWID_MARKER_CLOSE;
}