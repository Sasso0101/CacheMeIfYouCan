#include "mmio.h"
#include <graph.hpp>
#include <limits>

Reference::Reference(const CSR_local<uint32_t, float> *graph) : BFS_Impl(graph) {}

Reference::~Reference() {}

void Reference::BFS(vertex source, uint32_t *distances) {
  std::fill_n(distances, graph->nrows,
            std::numeric_limits<uint32_t>::max());
  std::vector<vertex> this_frontier = {};
  distances[source] = 0;
  this_frontier.push_back(source);
  #ifdef FRONTIER_DEBUG
  int i = 0;
  #endif
  while (!this_frontier.empty()) {
    #ifdef FRONTIER_DEBUG
    std::printf("frontier: %d size: %zu\n", i++, this_frontier.size());
    #endif
    std::vector<vertex> next_frontier;
    for (const auto &src : this_frontier) {
      for (uint64_t i = graph->row_ptr[src]; i < graph->row_ptr[src + 1]; i++) {
        vertex dst = graph->col_idx[i];
        if (distances[src] + 1 < distances[dst]) {
          distances[dst] = distances[src] + 1;
          next_frontier.push_back(dst);
        }
      }
    }
    std::swap(this_frontier, next_frontier);
  }
}

bool Reference::check_result(vertex source, uint32_t *distances) {
  return BFS_Impl::check_distances(source, distances);
}