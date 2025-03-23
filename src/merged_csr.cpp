#include "graph.hpp"
#include <limits>

#define DEGREE(vertex) merged[vertex]
#define DISTANCE(vertex) merged[vertex + 1]

MergedCSR::MergedCSR(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
    : BaseGraph(rowptr, col, N, M) {
  create_merged_csr(rowptr, col, N, M);
}

// Create merged CSR from CSR
void MergedCSR::create_merged_csr(eidType *rowptr, vidType *col, uint64_t N,
                                  uint64_t M) {
  merged = new eidType[M + 2 * N];
  eidType merged_index = 0;

#pragma omp parallel for schedule(static)
  for (vidType i = 0; i < N; i++) {
    eidType start = rowptr[i];
    // Add degree to start of neighbor list
    merged[merged_index++] = rowptr[i + 1] - rowptr[i];
    // Initialize distance
    merged[merged_index++] = std::numeric_limits<weight_type>::max();
    // Copy neighbors
    for (eidType j = start; j < rowptr[i + 1]; j++, merged_index++) {
      merged[merged_index] = rowptr[col[j]] + 2 * col[j];
    }
  }
  // Fix rowptr indices caused by adding the degree to the start of each
  // neighbor list
  for (vidType i = 0; i <= N; i++) {
    rowptr[i] = rowptr[i] + 2 * i;
  }
}

// Extract distances from merged CSR
void MergedCSR::compute_distances(weight_type *distances,
                                  vidType source) const {
#pragma omp parallel for simd schedule(static)
  for (vidType i = 0; i < N; i++) {
    distances[i] = DISTANCE(rowptr[i]);
    // Reset distance for next BFS
    DISTANCE(rowptr[i] + 1) = std::numeric_limits<weight_type>::max();
  }
  distances[source] = 0;
}

inline void MergedCSR::add_to_frontier(frontier &frontier, eidType v) const {
  frontier.push_back(v);
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
      eidType neighbor = merged[i];
      // If neighbor is not visited, add to frontier
      if (DISTANCE(neighbor) != std::numeric_limits<weight_type>::max()) {
        if (DEGREE(neighbor) != 1) {
          add_to_frontier(next_frontier, neighbor);
        }
        DISTANCE(neighbor) = distance;
      }
    }
  }
}

void MergedCSR::BFS(vidType source, weight_type *distances) {
  frontier this_frontier;
  eidType start = rowptr[source];

  add_to_frontier(this_frontier, start);
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