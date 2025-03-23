#include "graph.hpp"
#include "input.hpp"
#include <omp.h>
#include <string>

#define USAGE                                                                  \
  "Usage: %s <schema> [source] "                                               \
  "[heuristic|merged_csr|merged_csr_parents|bitmap|classic|reference] "        \
  "[check]\n The second argument specifies the source vertex. By default "     \
  "this is '0'. The third argument forces the usage of a specific algorithm. " \
  "By default a heuristic is used. The forth argument if the output of the "   \
  "algorithm should be verified against the reference implementation.\n"

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M, Algorithm algorithm) {
  switch (algorithm) {
  case Algorithm::MERGED_CSR_PARENTS:
    return new MergedCSR_Parents(rowptr, col, N, M);
    break;
  case Algorithm::MERGED_CSR:
    return new MergedCSR(rowptr, col, N, M);
    break;
  case Algorithm::BITMAP:
    return new Bitmap(rowptr, col, N, M);
    break;
  case Algorithm::CLASSIC:
    return new Classic(rowptr, col, N, M);
    break;
  case Algorithm::REFERENCE:
    return new Reference(rowptr, col, N, M);
    break;
  case Algorithm::HEURISTIC:
    if ((float)M / N < 10) { // Graph diameter heuristic
      return new MergedCSR(rowptr, col, N, M);
    } else {
      return new Bitmap(rowptr, col, N, M);
    }
    break;
  }
}

int main(const int argc, char **argv) {
  if (argc < 2 || argc > 6) {
    printf(USAGE, argv[0]);
    return 1;
  }
  std::string source = (argc >= 3) ? argv[2] : "0";
  std::string algorithm = (argc >= 4) ? argv[3] : "default";
  std::string problem = (argc >= 5) ? argv[4] : "distances";
  bool check = (argc == 6) ? true : false;

  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, algorithm, problem);
#pragma omp parallel
  {
#pragma omp master
    { printf("Number of threads: %d\n", omp_get_num_threads()); }
  }
  double t_start = omp_get_wtime();
  p.run(check, std::stoi(source));
  double t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);
  t_start = omp_get_wtime();
  p.run(check, 50);
  t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);
}