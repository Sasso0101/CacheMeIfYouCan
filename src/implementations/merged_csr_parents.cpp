#include <graph.hpp>
#include <omp.h>

#define VERTEX_ID(vertex) merged[vertex]
#define PARENT_ID(vertex) merged[vertex + 1]
#define DEGREE(vertex) merged[vertex + 2]

MergedCSR_Parents::MergedCSR_Parents(Graph *graph)
    : BFS_Impl(graph) {
  create_merged_csr();
}

MergedCSR_Parents::~MergedCSR_Parents() { delete[] merged; }

// Create merged CSR from CSR
void MergedCSR_Parents::create_merged_csr() {
  merged = new eidType[graph->M + 3 * graph->N];
  vidType merged_index = 0;
  for (vidType i = 0; i < graph->N; i++) {
    vidType start = graph->rowptr[i];
    // Add vertex ID to start of neighbor list
    merged[merged_index++] = i;
    // Add parent ID to start of neighbor list (initialized to -1)
    merged[merged_index++] = -1;
    // Add degree to start of neighbor list
    merged[merged_index++] = graph->rowptr[i + 1] - graph->rowptr[i];
    // Copy neighbors
    for (vidType j = start; j < graph->rowptr[i + 1]; j++, merged_index++) {
      merged[merged_index] = graph->rowptr[graph->col[j]] + 3 * graph->col[j];
    }
  }
  // Fix rowptr indices by adding offset caused by adding the degree to the
  // start of each neighbor list
  for (vidType i = 0; i <= graph->N; i++) {
    graph->rowptr[i] = graph->rowptr[i] + 3 * i;
  }
}
inline void MergedCSR_Parents::add_to_frontier(frontier &frontier,
                                               eidType v) const {
  frontier.push_back(v);
}

void MergedCSR_Parents::compute_parents(weight_type *parents,
                                        vidType source) const {
#pragma omp parallel for simd schedule(static)
  for (vidType i = 0; i < graph->N; i++) {
    parents[i] = merged[graph->rowptr[i] + 1];
  }
}

#pragma omp declare reduction(vec_add : std::vector<eidType> : omp_out.insert( \
        omp_out.end(), omp_in.begin(), omp_in.end()))

void MergedCSR_Parents::top_down_step(const frontier &this_frontier,
                                      frontier &next_frontier) {
#pragma omp parallel for reduction(vec_add : next_frontier)                    \
    schedule(static) if (this_frontier.size() > 50)
  for (const auto &v : this_frontier) {
    vidType end = v + merged[v + 2] + 3;
    for (eidType i = v + 3; i < end; i++) {
      eidType neighbor = merged[i];
      if (PARENT_ID(neighbor) == -1) {
        if (DEGREE(neighbor) != 1) {
          add_to_frontier(next_frontier, neighbor);
        }
        PARENT_ID(neighbor) = VERTEX_ID(v);
      }
    }
  }
}

void MergedCSR_Parents::BFS(vidType source, weight_type *parents) {
  frontier this_frontier;
  eidType start = graph->rowptr[source];
  eidType unexplored_edges = graph->M;
  Direction dir = Direction::TOP_DOWN;

  add_to_frontier(this_frontier, start);
  PARENT_ID(start) = source;
  while (!this_frontier.empty()) {
    frontier next_frontier;
    next_frontier.reserve(this_frontier.size());
    top_down_step(this_frontier, next_frontier);
    this_frontier = std::move(next_frontier);
  }
  compute_parents(parents, source);
}