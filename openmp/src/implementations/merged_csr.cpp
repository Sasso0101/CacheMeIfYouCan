#include "graph.hpp"
#include <limits>

#define DEGREE(vertex) merged_csr[vertex]
#define DISTANCE(vertex) merged_csr[vertex + 1]

MergedCSR_Distances::MergedCSR_Distances(const CSR_local<uint32_t, float>* graph) : BFS_Impl(graph) {
  create_merged_csr();
}

MergedCSR_Distances::~MergedCSR_Distances() { 
  delete[] merged_csr; 
  delete[] merged_rowptr;
}

// Create merged CSR from CSR
void MergedCSR_Distances::create_merged_csr() {
  merged_csr = new edge[graph->nnz + 2 * graph->nrows];
  merged_rowptr = new edge[graph->nrows + 1];
  edge merged_index = 0;

  for (vertex i = 0; i < graph->nrows; i++) {
    edge start = graph->row_ptr[i];
    // Add degree to start of neighbor list
    merged_csr[merged_index++] = graph->row_ptr[i + 1] - graph->row_ptr[i];
    // Initialize distance
    merged_csr[merged_index++] = std::numeric_limits<uint32_t>::max();
    // Copy neighbors
    for (edge j = start; j < graph->row_ptr[i + 1]; j++) {
      merged_csr[merged_index++] = graph->row_ptr[graph->col_idx[j]] + 2 * graph->col_idx[j];
    }
  }
  // Fix rowptr indices caused by adding the degree to the start of each
  // neighbor list
  for (vertex i = 0; i <= graph->nrows; i++) {
    merged_rowptr[i] = graph->row_ptr[i] + 2 * i;
  }
}

// Extract distances from merged CSR
void MergedCSR_Distances::compute_distances(edge *distances) const {
#pragma omp parallel for simd schedule(static)
  for (vertex i = 0; i < graph->nrows; i++) {
    distances[i] = DISTANCE(merged_rowptr[i]);
    // Reset distance for next BFS
    DISTANCE(merged_rowptr[i]) = std::numeric_limits<uint32_t>::max();
  }
}

#pragma omp declare reduction(vec_add                                          \
:frontier : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

void MergedCSR_Distances::top_down_step(const frontier &this_frontier,
                              frontier &next_frontier,
                              const uint32_t &distance) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    schedule(static) if (this_frontier.size() > 50)
  for (const auto &v : this_frontier) {
    edge end = v + 2 + DEGREE(v);
// Iterate over neighbors
    for (edge i = v + 2; i < end; i++) {
      edge neighbor = merged_csr[i];
      // If neighbor is not visited, add to frontier
      if (DISTANCE(neighbor) == std::numeric_limits<uint32_t>::max()) {
        if (DEGREE(neighbor) != 1) {
          next_frontier.push_back(neighbor);
        }
        DISTANCE(neighbor) = distance;
      }
    }
  }
}

void MergedCSR_Distances::BFS(vertex source, uint32_t *distances) {
  frontier this_frontier;
  edge start = merged_rowptr[source];

  this_frontier.push_back(start);
  DISTANCE(start) = 0;
  uint32_t distance = 1;
  while (!this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    top_down_step(this_frontier, next_frontier, distance);
    distance++;
    this_frontier = std::move(next_frontier);
  }
  compute_distances(distances);
}

bool MergedCSR_Distances::check_result(vertex source, uint32_t *distances) {
  return BFS_Impl::check_distances(source, distances);
}