#pragma once
#include <cstdint>
#include <vector>

// #define DBG_THREAD_BALANCE
#ifdef DBG_THREAD_BALANCE
  #include <omp.h>
  #include <cstdio>
  #include <utility>
#endif

typedef uint32_t vidType;
typedef uint64_t eidType;
typedef uint32_t weight_type;

enum class Direction { TOP_DOWN, BOTTOM_UP };
typedef std::vector<vidType> frontier;

typedef struct {
  vidType parent;
  vidType node;
} pruned;

class BaseGraph {
public:
  #ifdef DBG_THREAD_BALANCE
    std::vector<uint32_t> thread_balance_niter;
    std::vector<uint32_t> thread_balance_nwrites;
    std::vector<std::pair<weight_type, Direction>> td_bu_switches;
    BaseGraph() {
      int nthreads = omp_get_max_threads();
      thread_balance_niter.resize(nthreads, 0);
      thread_balance_nwrites.resize(nthreads, 0);
    }
  #else
    BaseGraph() {}
  #endif
  virtual void BFS(vidType source, weight_type *distances) = 0;
};