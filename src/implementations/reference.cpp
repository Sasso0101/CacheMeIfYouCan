#include <implementation.hpp>

Reference::Reference(Graph *graph) : BFS_Impl(graph) {}

Reference::~Reference() {}

void Reference::BFS(vidType source, weight_type *distances) {
  std::vector<vidType> this_frontier = {};
  distances[source] = 0;
  this_frontier.push_back(source);
  while (!this_frontier.empty()) {
    std::vector<vidType> next_frontier;
    for (const auto &src : this_frontier) {
      for (uint64_t i = graph->row_ptr[src]; i < graph->row_ptr[src + 1]; i++) {
        vidType dst = graph->col_idx[i];
        if (distances[src] + 1 < distances[dst]) {
          distances[dst] = distances[src] + 1;
          next_frontier.push_back(dst);
        }
      }
    }
    std::swap(this_frontier, next_frontier);
  }
}

bool Reference::check_result(vidType source, weight_type *distances) {
  return BFS_Impl::check_distances(source, distances);
}