#include "input.hpp"
#include <profiling.hpp>

int main(const int argc, char **argv) {
  std::ifstream in("schema/" + std::string(argv[1]));
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, complete::initialize_graph);
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