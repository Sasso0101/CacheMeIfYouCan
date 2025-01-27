#include "input.hpp"
#include <omp.h>
#include <profiling.hpp>
#include <complete.hpp>

#define USAGE "Usage: %s <schema> [small|large|very_large]\nThe second argument forces the usage of a specific algorithm. By default this is choosen dynamically.\n"

int main(const int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    printf(USAGE, argv[0]);
    return 1;
  }
  std::string algorithm = (argc == 3) ? argv[2] : "default";

  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, complete::initialize_graph, algorithm);
  #pragma omp parallel
  {
    #pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
  }
  double t_start = omp_get_wtime();
  p.check();
  double t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);
}