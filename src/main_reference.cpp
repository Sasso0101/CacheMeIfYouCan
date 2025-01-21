#include "input.hpp"
#include <omp.h>
#include <profiling.hpp>
#include <complete.hpp>

int main(const int argc, char **argv) {
  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, reference::initialize_graph, std::string(""));
  LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    #pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
    LIKWID_MARKER_THREADINIT;
  }
  #pragma omp parallel
  {
    LIKWID_MARKER_START("bfs");
  }
  double t_start = omp_get_wtime();
  p.run();
  double t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);
  #pragma omp parallel
  {
    LIKWID_MARKER_STOP("bfs");
  }
  LIKWID_MARKER_CLOSE;
}