#include <graph.hpp>
#include <omp.h>

#define VERTEX_ID(vertex) merged_csr[vertex]
#define PARENT_ID(vertex) merged_csr[vertex + 1]
#define DEGREE(vertex) merged_csr[vertex + 2]

MergedCSR_Parents::MergedCSR_Parents(const CSR_local<uint32_t, float>* graph) : BFS_Impl(graph) {
  create_merged_csr();
}

MergedCSR_Parents::~MergedCSR_Parents() { delete[] merged_csr; }

// Create merged CSR from CSR
void MergedCSR_Parents::create_merged_csr() {
  merged_csr = new edge[graph->nnz + 3 * graph->nrows];
  merged_rowptr = new edge[graph->nnz];

  vertex merged_index = 0;
  for (vertex i = 0; i < graph->nrows; i++) {
    vertex start = graph->row_ptr[i];
    // Add vertex ID to start of neighbor list
    merged_csr[merged_index++] = i;
    // Add parent ID to start of neighbor list (initialized to -1)
    merged_csr[merged_index++] = -1;
    // Add degree to start of neighbor list
    merged_csr[merged_index++] = graph->row_ptr[i + 1] - graph->row_ptr[i];
    // Copy neighbors
    for (vertex j = start; j < graph->row_ptr[i + 1]; j++) {
      merged_csr[merged_index++] =
          graph->row_ptr[graph->col_idx[j]] + 3 * graph->col_idx[j];
    }
  }
  // Fix rowptr indices by adding offset caused by adding the degree to the
  // start of each neighbor list
  for (vertex i = 0; i <= graph->nrows; i++) {
    merged_rowptr[i] = graph->row_ptr[i] + 3 * i;
  }
}

void MergedCSR_Parents::compute_parents(uint32_t *parents) const {
#pragma omp parallel for simd schedule(static)
  for (vertex i = 0; i < graph->nrows; i++) {
    parents[i] = merged_csr[merged_rowptr[i] + 1];
  }
}

#pragma omp declare reduction(vec_add : std::vector<edge> : omp_out.insert( \
        omp_out.end(), omp_in.begin(), omp_in.end()))

void MergedCSR_Parents::top_down_step(const frontier &this_frontier,
                                      frontier &next_frontier) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    schedule(static) if (this_frontier.size() > 50)
  for (const auto &v : this_frontier) {
    vertex end = v + DEGREE(v) + 3;
    for (edge i = v + 3; i < end; i++) {
      edge neighbor = merged_csr[i];
      if (PARENT_ID(neighbor) == -1.0) {
        if (DEGREE(neighbor) != 1) {
          next_frontier.push_back(neighbor);
        }
        PARENT_ID(neighbor) = VERTEX_ID(v);
      }
    }
  }
}

void MergedCSR_Parents::BFS(vertex source, uint32_t *parents) {
  frontier this_frontier = {};
  edge start = merged_rowptr[source];

  this_frontier.push_back(start);
  PARENT_ID(start) = source;
  while (!this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    top_down_step(this_frontier, next_frontier);
    this_frontier = std::move(next_frontier);
  }
  compute_parents(parents);
}

bool MergedCSR_Parents::check_result(vertex source, uint32_t *parents) {
  return BFS_Impl::check_parents(source, parents);
}