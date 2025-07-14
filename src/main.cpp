#include "../distributed_mmio/include/mmio.h"
#include "implementation.hpp"

#include <algorithm>
#include <limits>
#include <omp.h>
#include <random>
#include <string>

#define USAGE                                                                  \
  "Usage: %s <matrix> <runs> <implementation> <check>\n"                       \
  "Runs BFS implementations. \n\n"                                             \
  "Mandatory arguments:\n"                                                     \
  "<matrix>\t : path to matrix\n"                                              \
  "<runs>\t\t : integer. Number of runs to execute (default: 1) \n"            \
  "<algorithm>\t : 'merged_csr_parents', 'merged_csr', 'bitmap', 'classic', "  \
  "'reference', 'heuristic' (default: 'heuristic') \n"                         \
  "<check>\t\t : 'true', false'. Checks correctness of the result ('false' "   \
  "by "                                                                        \
  "default)\n"

// Seed used for picking source vertices
// Using same seed as in GAP benchmark for reproducible experiments
// https://github.com/sbeamer/gapbs/blob/b5e3e19c2845f22fb338f4a4bc4b1ccee861d026/src/util.h#L22
#define SEED 27491095

BFS_Impl *initialize_BFS(const char *filename, std::string algorithm) {
  Graph *csr_matrix =
      Distr_MMIO_CSR_local_read<uint32_t, float>(filename, false, NULL);

  if (algorithm == "merged_csr_parents") {
    return new MergedCSR_Parents(csr_matrix);
  } else if (algorithm == "merged_csr") {
    return new MergedCSR(csr_matrix);
  } else if (algorithm == "bitmap") {
    return new Bitmap(csr_matrix);
  } else if (algorithm == "classic") {
    return new Classic(csr_matrix);
  } else if (algorithm == "reference") {
    return new Reference(csr_matrix);
  } else {
    if ((float)(csr_matrix->nnz) / csr_matrix->nrows <
        10) { // Graph diameter heuristic
      return new MergedCSR(csr_matrix);
    } else {
      return new Bitmap(csr_matrix);
    }
  }
}

int main(const int argc, char **argv) {
  if (argc < 2 || argc > 5) {
    printf(USAGE, argv[0]);
    return 1;
  }
  vidType num_runs = 1;
  bool check = false;

  if (argc > 2) {
    num_runs = std::stoi(argv[2]);
  }
  std::string algorithm = "heuristic";
  if (argc > 3) {
    algorithm = argv[3];
  }
  if (argc > 4) {
    std::string check_str = argv[4];
    if (check_str == "true") {
      check = true;
    }
  }

#pragma omp parallel
  {
#pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
  }

  BFS_Impl *bfs = initialize_BFS(argv[1], algorithm);

  weight_type *result = new weight_type[bfs->graph->nrows];

  // Random number generator for source picking
  std::mt19937_64 rng(SEED); // Seed for reproducibility
  std::uniform_int_distribution<vidType> udist(0, bfs->graph->nrows - 1);
  double t_start, t_end;

  for (int run = 0; run < num_runs; ++run) {
    // Pick a random source vertex
    vidType random_source;
    do {
      random_source = udist(rng);
    } while (bfs->graph->row_ptr[random_source] ==
             bfs->graph->row_ptr[random_source + 1]);

    // Initialize result vector
    std::fill_n(result, bfs->graph->nrows,
                std::numeric_limits<weight_type>::max());

    t_start = omp_get_wtime();
    bfs->BFS(random_source, result);
    t_end = omp_get_wtime();

    printf("run=%d,%f\n", run + 1, t_end - t_start);

    if (check) {
      bfs->check_result(random_source, result);
    }
  }
  delete[] result;
  return 0;
}