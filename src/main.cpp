#include "input.hpp"
#include <omp.h>
// #include <profiling.hpp>
#include <complete.hpp>
#include <string>

#define USAGE "Usage: %s <schema> [source] [small|large|classic] [distances|parents] [check]\n The second argument specifies the source vertex. By default this is '0'. The third argument forces the usage of a specific algorithm. By default this is choosen dynamically. The forth argument specifies the output of the algorithm. By default this is 'distances'.\n"

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
  // LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
    // int thread_id = omp_get_thread_num();
    #pragma omp master
    {
      printf("Number of threads: %d\n", omp_get_num_threads());
    }
    // LIKWID_MARKER_THREADINIT;
  }
  // #pragma omp parallel
  // {
  //   LIKWID_MARKER_START("bfs");
  // }

  #ifdef DBG_CACHE
    int perf_ctl_fd;
    int perf_ctl_ack_fd;
    char ack[5];
    perf_ctl_fd = atoi(getenv("PERF_CTL_FD"));
    perf_ctl_ack_fd = atoi(getenv("PERF_CTL_ACK_FD"));

    write(perf_ctl_fd, "enable\n", 8);
    read(perf_ctl_ack_fd, ack, 5);
    assert(strcmp(ack, "ack\n") == 0);
  #endif
  double t_start = omp_get_wtime();
  p.run(check, std::stoi(source));
  double t_end = omp_get_wtime();
  printf("Runtime: %f\n", t_end - t_start);

  #ifdef DBG_CACHE
    write(perf_ctl_fd, "disable\n", 9);
    read(perf_ctl_ack_fd, ack, 5);
    assert(strcmp(ack, "ack\n") == 0);
  #endif

  // #pragma omp parallel
  // {
  //   LIKWID_MARKER_STOP("bfs");
  // }
  // LIKWID_MARKER_CLOSE;
}