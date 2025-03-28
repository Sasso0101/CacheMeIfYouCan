#include "graph.hpp"
#include <CLI/CLI.hpp>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <random>
#include <string>

BFS_Impl *initialize_BFS(std::string filename, Implementation impl) {
  std::string path = "schemas/" + filename;
  Graph *graph = new Graph(path);

  switch (impl) {
  case merged_csr_parents:
    return new MergedCSR_Parents(graph);
  case merged_csr:
    return new MergedCSR(graph);
  case merged_csr_1:
    return new MergedCSR_1(graph);
  case bitmap:
    return new Bitmap(graph);
  case classic:
    return new Classic(graph);
  case reference:
    return new Reference(graph);
  case heuristic:
    if ((float)(graph->M) / graph->N < 10) { // Graph diameter heuristic
      return new MergedCSR(graph);
    } else {
      return new Bitmap(graph);
    }
  }
  fprintf(stderr, "Error: Unsupported implementation.\n");
  exit(EXIT_FAILURE);
}

std::vector<vidType> pick_random_vertices(Graph *graph, int count, int seed) {
  std::vector<vidType> vertices;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(0, graph->N - 1);
  for (int i = 0; i < count; i++) {
    vertices.push_back(dis(gen));
  }
  return vertices;
}

std::map<std::string, Implementation> map {
    {"merged_csr_parents", merged_csr_parents},
    {"merged_csr", merged_csr},
    {"merged_csr_1", merged_csr_1},
    {"bitmap", bitmap},
    {"classic", classic},
    {"reference", reference},
    {"heuristic", heuristic}
};

double run_bfs(BFS_Impl *bfs, vidType source, bool check) {
  weight_type *result = new weight_type[bfs->graph->N];
  std::fill_n(result, bfs->graph->N, std::numeric_limits<weight_type>::max());
  double t_start = omp_get_wtime();
  bfs->BFS(source, result);
  double t_end = omp_get_wtime();
  if (check) {
    bfs->check_result(source, result);
  }
  delete [] result;
  return t_end - t_start;
}

int main(const int argc, char **argv) {
  CLI::App app{"Optimized BFS"};
  argv = app.ensure_utf8(argv);

  std::string schema;
  Implementation impl = heuristic;
  int seed = 0;
  int runs = 1;
  bool check = false;

  app.add_option("-f,--schema", schema, "Path to JSON schema of dataset")->required();
  app.add_option("-i,--implementation", impl, "Implementation to use: 'merged_csr_parents', 'merged_csr', 'bitmap', 'classic', 'reference', 'heuristic' ('heuristic' by default)")->transform(CLI::CheckedTransformer(map, CLI::ignore_case));

  app.add_option("-s,--seed", seed, "Seed for random number generator");
  app.add_option("-r,--runs", runs, "Number of runs");
  app.add_flag("-c,--check", check, "Check correctness of the result");

  CLI11_PARSE(app, argc, argv);
#pragma omp parallel
  {
#pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
  }
  double t_start = omp_get_wtime();
  BFS_Impl *bfs = initialize_BFS(schema, impl);
  double t_end = omp_get_wtime();

  printf("Algorithm: %s\n", schema.c_str());
  printf("Initialization: %f\n", t_end - t_start);

  std::vector<vidType> sources = pick_random_vertices(bfs->graph, runs, seed);
  run_bfs(bfs, sources[0], false); // Dry run to warm up
  printf("Runtime: ");
  for (int i = 0; i < runs; i++) {
    double time = run_bfs(bfs, sources[i], check);
    printf("%f ", time);
  }
  printf("\n");
}