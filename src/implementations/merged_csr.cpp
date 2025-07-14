#include "implementation.hpp"
#include <limits>

#define DEGREE(vertex) merged_csr[vertex]
#define DISTANCE(vertex) merged_csr[vertex + 1]

MergedCSR::MergedCSR(Graph *graph) : BFS_Impl(graph) {
  create_merged_csr();
}

MergedCSR::~MergedCSR() { delete[] merged_csr; }

// Create merged CSR from CSR
void MergedCSR::create_merged_csr() {
  merged_csr = new eidType[graph->nnz + 2 * graph->nrows];
  merged_rowptr = new eidType[graph->nrows];
  eidType merged_index = 0;

  for (vidType i = 0; i < graph->nrows; i++) {
    eidType start = graph->row_ptr[i];
    // Add degree to start of neighbor list
    merged_csr[merged_index++] = graph->row_ptr[i + 1] - graph->row_ptr[i];
    // Initialize distance
    merged_csr[merged_index++] = std::numeric_limits<weight_type>::max();
    // Copy neighbors
    for (eidType j = start; j < graph->row_ptr[i + 1]; j++) {
      merged_csr[merged_index++] = graph->row_ptr[graph->col_idx[j]] + 2 * graph->col_idx[j];
    }
  }
  // Fix rowptr indices caused by adding the degree to the start of each
  // neighbor list
  for (vidType i = 0; i <= graph->nrows; i++) {
    merged_rowptr[i] = graph->row_ptr[i] + 2 * i;
  }
}

// Extract distances from merged CSR
void MergedCSR::compute_distances(weight_type *distances,
                                  vidType source) const {
#pragma omp parallel for simd schedule(static)
  for (vidType i = 0; i < graph->nrows; i++) {
    distances[i] = DISTANCE(merged_rowptr[i]);
    // Reset distance for next BFS
    DISTANCE(merged_rowptr[i] + 1) = std::numeric_limits<weight_type>::max();
  }
  distances[source] = 0;
}

#pragma omp declare reduction(vec_add                                          \
:frontier : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

void MergedCSR::top_down_step(const frontier &this_frontier,
                              frontier &next_frontier,
                              const weight_type &distance) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    schedule(static) if (this_frontier.size() > 50)
  for (const auto &v : this_frontier) {
    eidType end = v + 2 + DEGREE(v);
// Iterate over neighbors
#pragma omp simd
    for (eidType i = v + 2; i < end; i++) {
      eidType neighbor = merged_csr[i];
      // If neighbor is not visited, add to frontier
      if (DISTANCE(neighbor) == std::numeric_limits<weight_type>::max()) {
        if (DEGREE(neighbor) != 1) {
          next_frontier.push_back(neighbor);
        }
        DISTANCE(neighbor) = distance;
      }
    }
  }
}

void MergedCSR::BFS(vidType source, weight_type *distances) {
  frontier this_frontier;
  eidType start = merged_rowptr[source];

  this_frontier.push_back(start);
  DISTANCE(start) = 0;
  weight_type distance = 1;
  while (!this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    top_down_step(this_frontier, next_frontier, distance);
    distance++;
    this_frontier = std::move(next_frontier);
  }
  compute_distances(distances, source);
}

bool MergedCSR::check_result(vidType source, weight_type *distances) {
  return BFS_Impl::check_distances(source, distances);
}