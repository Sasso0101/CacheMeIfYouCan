#include "input.hpp"
#include <omp.h>
#include <profiling.hpp>
#include <complete.hpp>

#define USAGE "Usage: %s <schema> [small|large|very_large] [distances|parents]\nThe second argument forces the usage of a specific algorithm. By default this is choosen dynamically. The third argument specifies the output of the algorithm.\n"

int main(const int argc, char **argv) {
  if (argc < 2 || argc > 4) {
    printf(USAGE, argv[0]);
    return 1;
  }
  std::string algorithm = (argc >= 3) ? argv[2] : "default";
  std::string output = (argc == 4) ? argv[3] : "distances";

  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, algorithm, output);
  LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
    // int thread_id = omp_get_thread_num();
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
  p.run(true, 0);
  double t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);
  #pragma omp parallel
  {
    LIKWID_MARKER_STOP("bfs");
  }
  LIKWID_MARKER_CLOSE;
}