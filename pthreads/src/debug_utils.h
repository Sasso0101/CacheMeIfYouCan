#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

#include "config.h"
#include "frontier.h"
#include "merged_csr.h"
#include "mmio_c_wrapper.h"
#include <stdint.h>
#include <stdio.h>

void print_chunk_counts(const Frontier *f) {
  printf("Chunk counts: ");
  for (int i = 0; i < MAX_THREADS; i++) {
    printf("%5d", f->thread_chunk_counts[i]);
  }
  printf("\n");
}

void print_sources(const mmio_csr_u32_f32_t *graph, const uint32_t *sources, int runs) {
  printf("Sources: ");
  for (int i = 0; i < runs; i++) {
    printf("%u (%u)\n", sources[i], graph->row_ptr[sources[i]]);
  }
  printf("\n");
}

void print_graph(const mmio_csr_u32_f32_t *graph) {
  for (uint32_t i = 0; i < graph->nrows; i++) {
    printf("Vertex %u: ", i);
    for (uint32_t j = graph->row_ptr[i]; j < graph->row_ptr[i + 1]; j++) {
      printf("%u ", graph->col_idx[j]);
    }
    printf("\n");
  }
}

void print_merged_csr(const MergedCSR *merged_csr) {
  for (uint32_t i = 0; i < merged_csr->num_vertices; i++) {
    printf("Vertex %u (%lu): ", i, (uint64_t)merged_csr->row_ptr[i]);
    for (mer_t j = merged_csr->row_ptr[i]; j < merged_csr->row_ptr[i + 1];
         j++) {
      printf("%lu ", (uint64_t)merged_csr->merged[j]);
    }
    printf("\n");
  }
}

void print_distances(const mmio_csr_u32_f32_t *graph, const uint32_t *distances) {
  printf("Distances:\n");
  for (uint32_t i = 0; i < graph->nrows; i++) {
    printf("Vertex %u: %u\n", i, distances[i]);
  }
}

void print_time(const double elapsed[], int runs) {
  double total = 0;
  for (int i = 0; i < runs; i++) {
    total += elapsed[i];
  }
  printf("Average time: %14.5f\n", total / runs);
}

int check_bfs_correctness(const mmio_csr_u32_f32_t *graph, const uint32_t *distances,
                          uint32_t source) {
  uint32_t n = graph->nrows;

  // Basic sanity check for source index
  if (source >= n) {
    printf("Error: Source vertex %u is out of bounds (num_vertices=%u)\n",
           source, n);
    return 0;
  }

  // --- Check 1: Source distance must be 0 ---
  if (distances[source] != 0) {
    printf("Error: Distance to source vertex %u must be 0 (got %u)\n", source,
           distances[source]);
    return 0;
  }

  // --- Iterate through all vertices to check properties ---
  for (uint32_t u = 0; u < n; u++) {
    uint32_t dist_u = distances[u];
    uint32_t start = graph->row_ptr[u];
    uint32_t end = graph->row_ptr[u + 1];
    bool found_predecessor = false; // Only relevant if dist_u > 0

    // --- Check 2: Unreachable nodes ---
    if (dist_u == UINT32_MAX) {
      // If a node is marked unreachable, none of its neighbors should be
      // reachable. (Because if a neighbor v was reachable, u would be reachable
      // too in an undirected graph)
      for (uint32_t i = start; i < end; i++) {
        uint32_t v = graph->col_idx[i];
        if (v >= n) { // Basic bounds check for neighbor index
          printf("Error: Neighbor index %u for vertex %u is out of bounds.\n",
                 v, u);
          return 0;
        }
        if (distances[v] != UINT32_MAX) {
          printf("Error: Vertex %u is marked unreachable (dist=%u), but has a "
                 "reachable neighbor %u (dist=%u)\n",
                 u, dist_u, v, distances[v]);
          return 0;
        }
      }
      // If u is the source AND unreachable, Check 1 already failed.
    }
    // --- Check 3: Reachable nodes (including source) ---
    else {
      // Check 3a: Nodes other than source cannot have distance 0
      if (dist_u == 0 && u != source) {
        printf(
            "Error: Vertex %u has distance 0 but is not the source vertex %u\n",
            u, source);
        return 0;
      }

      // Check 3b: For reachable nodes, check neighbors
      for (uint32_t i = start; i < end; i++) {
        uint32_t v = graph->col_idx[i];
        if (v >= n) { // Basic bounds check for neighbor index
          printf("Error: Neighbor index %u for vertex %u is out of bounds.\n",
                 v, u);
          return 0;
        }
        uint32_t dist_v = distances[v];

        // Rule 1 (Edge Relaxation): dist[v] must be <= dist[u] + 1
        // This handles the transition from reachable to unreachable (UINT32_MAX
        // is large) and checks for jumps that are too large.
        if (dist_v > dist_u + 1) {
          printf("Error: Edge {%u, %u} violates distance rule. dist[%u]=%u, "
                 "dist[%u]=%u. Neighbor distance is too large.\n",
                 u, v, u, dist_u, v, dist_v);
          return 0;
        }

        // Rule 2 (Predecessor Check preparation): If u is not the source,
        // track if we find a neighbor v such that dist[v] == dist[u] - 1.
        if (dist_u > 0 && dist_v == dist_u - 1) {
          found_predecessor = true;
        }
      } // End neighbor loop

      // Check 3c: Predecessor Check verification
      // Every reachable node 'u' (except the source) must have at least one
      // neighbor 'v' with distance exactly dist[u] - 1.
      if (dist_u > 0 && !found_predecessor) {
        printf("Error: Reachable vertex %u (dist=%u) has no predecessor "
               "neighbor with distance %u.\n",
               u, dist_u, dist_u - 1);
        return 0;
      }
    } // End reachable node checks
  } // End main vertex loop

  // If all checks passed
  printf("BFS result verification passed for source vertex %u.\n", source);
  return 1;
}

#endif // DEBUG_UTILS_H