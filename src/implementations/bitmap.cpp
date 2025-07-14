#include "implementation.hpp"

#define IS_VISITED(i) (visited[i])

inline void Bitmap::add_to_frontier(bool *frontier, vidType v) {
  frontier[v] = true;
  visited[v] = true;
}

Bitmap::Bitmap(Graph *graph)
    : BFS_Impl(graph), this_frontier(new bool[graph->nrows]),
      next_frontier(new bool[graph->nrows]), visited(new bool[graph->nrows]) {
#pragma omp parallel for schedule(static)
  for (eidType i = 0; i < graph->nrows; i++) {
    this_frontier[i] = false;
    next_frontier[i] = false;
    visited[i] = false;
  }
}

Bitmap::~Bitmap() {
  delete[] this_frontier;
  delete[] next_frontier;
  delete[] visited;
}

void Bitmap::bottom_up_step(const bool *this_frontier, bool *next_frontier) {
#pragma omp parallel for schedule(static)
  for (vidType i = 0; i < graph->nrows; i++) {
    if (!IS_VISITED(i)) {
      for (eidType j = graph->row_ptr[i]; j < graph->row_ptr[i + 1]; j++) {
        vidType neighbor = graph->col_idx[j];
        if (this_frontier[neighbor] == true) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, i);
          break;
        }
      }
    }
  }
}

void Bitmap::top_down_step(const bool *this_frontier, bool *next_frontier) {
#pragma omp parallel for schedule(static)
  for (int v = 0; v < graph->nrows; v++) {
    if (this_frontier[v] == true) {
      eidType end = graph->row_ptr[v + 1];
#pragma omp simd
      for (eidType i = graph->row_ptr[v]; i < end; i++) {
        vidType neighbor = graph->col_idx[i];
        if (!IS_VISITED(neighbor)) {
          add_to_frontier(next_frontier, neighbor);
        }
      }
    }
  }
}

void Bitmap::BFS(vidType source, weight_type *distances) {
  eidType unexplored_edges = graph->nnz;
  eidType unvisited_vertices = graph->nrows;
  Direction dir = Direction::TOP_DOWN;
  add_to_frontier(this_frontier, source);
  eidType edges_frontier = graph->row_ptr[source + 1] - graph->row_ptr[source];
  vidType vertices_frontier = 1;
  distances[source] = 0;
  weight_type distance = 1;

  do {
    if (dir == Direction::BOTTOM_UP && vertices_frontier < graph->nrows / BETA) {
      dir = Direction::TOP_DOWN;
    } else if (dir == Direction::TOP_DOWN &&
               edges_frontier > unexplored_edges / ALPHA) {
      dir = Direction::BOTTOM_UP;
    }
    unexplored_edges -= edges_frontier;
    unvisited_vertices -= vertices_frontier;
    edges_frontier = 0;
    vertices_frontier = 0;
    if (dir == Direction::TOP_DOWN) {
      top_down_step(this_frontier, next_frontier);
    } else {
      bottom_up_step(this_frontier, next_frontier);
    }
#pragma omp parallel for reduction(+ : edges_frontier, vertices_frontier)      \
    schedule(static)
    for (vidType i = 0; i < graph->nrows; i++) {
      this_frontier[i] = false;
      if (next_frontier[i] == true) {
        edges_frontier += graph->row_ptr[i + 1] - graph->row_ptr[i];
        vertices_frontier += 1;
        distances[i] = distance;
      }
    }
    if (vertices_frontier == 0) {
      break;
    }
    std::swap(this_frontier, next_frontier);
    distance++;
  } while (true);
#pragma omp parallel for schedule(static)
  for (vidType i = 0; i < graph->nrows; i++) {
    this_frontier[i] = false;
    visited[i] = false;
  }
}

bool Bitmap::check_result(vidType source, weight_type *distances) {
  return BFS_Impl::check_distances(source, distances);
}