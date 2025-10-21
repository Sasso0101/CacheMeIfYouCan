#include "graph.hpp"
#include <omp.h>
#include <random>
#include <sys/types.h>

#ifdef USE_PAPI
extern "C" {
#include <papi.h>
}
void handle_error(int retval) {
  printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
  exit(1);
}
#endif

#define USAGE                                                                  \
  "Usage: %s <dataset> <runs> <implementation> <check> <source> \nRuns BFS "   \
  "implementations. \n\nMandatory arguments:\n  <dataset>\t path to dataset "  \
  "\n  <runs>\t\t : integer. Number of runs (1 by default) \n  <source>\t : "  \
  "integer. Source vertex ID (64 randomly generated vertices by default) \n "  \
  " <algorithm>\t : 'merged_csr_parents', 'merged_csr_distances', "            \
  "'reference' ('reference' by default) \n  <check>\t : 'true', false'. "      \
  "Checks correctness of the result ('false' by default)\n"

// A constant seed for the random number generator, equal to kRandSeed in GAPBS.
const int kRandSeed = 27491095;

// UniDist is a utility class from the GAP Benchmark Suite for generating
// uniformly distributed random numbers within a specific range.
template <typename NodeID_, typename rng_t_,
          typename uNodeID_ = typename std::make_unsigned<NodeID_>::type>
class UniDist {
public:
  UniDist(NodeID_ max_value, rng_t_ &rng) : rng_(rng) {
    uNodeID_ u_max_value = static_cast<uNodeID_>(max_value);
    no_mod_ = (rng_.max() == u_max_value);
    mod_ = u_max_value + 1;
    uNodeID_ remainder_sub_1 = rng_.max() % mod_;
    if (remainder_sub_1 == mod_ - 1)
      cutoff_ = 0;
    else
      cutoff_ = rng_.max() - remainder_sub_1;
  }

  NodeID_ operator()() {
    uNodeID_ rand_num = rng_();
    if (no_mod_)
      return rand_num;
    if (cutoff_ != 0) {
      while (rand_num >= cutoff_)
        rand_num = rng_();
    }
    return rand_num % mod_;
  }

private:
  rng_t_ &rng_;
  bool no_mod_;
  uNodeID_ mod_;
  uNodeID_ cutoff_;
};

// Generates a vector of random source vertices for a given graph.
// It ensures that selected vertices have an out-degree greater than zero.
void generate_random_sources(const CSR_local<uint32_t, float> *graph,
                        size_t num_sources, std::vector<uint32_t>& sources) {
  std::mt19937_64 rng(kRandSeed);
  UniDist<uint32_t, std::mt19937_64> udist(graph->nrows - 1, rng);

  while (sources.size() < num_sources) {
    uint32_t source = udist();
    // Ensure the source has outgoing edges
    if ((graph->row_ptr[source + 1] - graph->row_ptr[source]) > 0) {
      sources.push_back(source);
    }
  }
}

int main(const int argc, char **argv) {
  if (argc < 2 || argc > 6) {
    printf(USAGE, argv[0]);
    return 1;
  }
  std::vector<uint32_t> sources = {};
  bool check = false;
  int runs = 1;

  std::string algo_str = "reference";
  if (argc > 3) {
    algo_str = std::string(argv[3]);
  }

  double t_start = omp_get_wtime();
  CSR_local<uint32_t, float> *graph =
      Distr_MMIO_CSR_local_read<uint32_t, float>(argv[1], false);

  BFS_Impl *bfs;
  if (algo_str == "merged_csr_parents") {
    printf("Using Merged CSR with Parents implementation\n");
    bfs = new MergedCSR_Parents(graph);
  } else if (algo_str == "merged_csr_distances") {
    printf("Using Merged CSR with Distances implementation\n");
    bfs = new MergedCSR_Distances(graph);
  } else {
    printf("Using Reference implementation\n");
    bfs = new Reference(graph);
  }
  double t_end = omp_get_wtime();

  printf("Initialization: %f\n", t_end - t_start);

  if (argc > 2) {
    runs = std::stoi(argv[2]);
  }

  if (argc > 4) {
    check = true;
  }

  if (argc > 5) {
    sources.insert(sources.end(), runs, std::stoi(argv[5]));
  } else {
    generate_random_sources(bfs->graph, runs, sources);
  }

#pragma omp parallel
  {
#pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
  }

  uint32_t *result = new uint32_t[bfs->graph->nrows];

  for (uint32_t i = 0; i < sources.size(); i++) {
#ifndef USE_PAPI
    t_start = omp_get_wtime();
#endif
    #ifdef USE_PAPI
    if (i == 2) { // skip first two iterations
      #pragma omp parallel
      {
        int retval = PAPI_hl_region_begin("computation");
        if ( retval != PAPI_OK )
          handle_error(retval);
      }
    }
    #endif
    bfs->BFS(sources[i], result);
#ifndef USE_PAPI
    t_end = omp_get_wtime();
    printf("run_id=%d,threads=%d,source=%d,%.4f\n", i, omp_get_max_threads(),
           sources[i], t_end - t_start);
    if (check) {
      printf("Checking result for source %d\n", sources[i]);
      bfs->check_result(sources[i], result);
    }
#endif
  }
#ifdef USE_PAPI
  #pragma omp parallel
  {
    int retval = PAPI_hl_region_end("computation");
    if (retval != PAPI_OK)
      handle_error(retval);
  }
#endif
  delete[] result;
  delete graph;
  delete bfs;
  return 0;
}