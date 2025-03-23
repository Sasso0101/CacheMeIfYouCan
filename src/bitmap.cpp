#include "graph.hpp"

#define IS_VISITED(i) (visited[i])

inline void Bitmap::add_to_frontier(bool *frontier, vidType v) {
  frontier[v] = true;
  visited[v] = true;
}

Bitmap::Bitmap(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
    : BaseGraph(rowptr, col, N, M) {
  bool *this_frontier = new bool[N];
  bool *next_frontier = new bool[N];
  bool *visited = new bool[N];
#pragma omp parallel for schedule(static)
  for (eidType i = 0; i < N; i++) {
    this_frontier[i] = false;
    next_frontier[i] = false;
    visited[i] = false;
  }
}

void Bitmap::bottom_up_step(const bool *this_frontier, bool *next_frontier) {
#pragma omp parallel for schedule(static)
  for (vidType i = 0; i < N; i++) {
    if (!IS_VISITED(i)) {
      for (eidType j = rowptr[i]; j < rowptr[i + 1]; j++) {
        vidType neighbor = col[j];
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
  for (int v = 0; v < N; v++) {
    if (this_frontier[v] == true) {
      eidType end = rowptr[v + 1];
#pragma omp simd
      for (eidType i = rowptr[v]; i < end; i++) {
        vidType neighbor = col[i];
        if (!IS_VISITED(neighbor)) {
          add_to_frontier(next_frontier, neighbor);
        }
      }
    }
  }
}

void Bitmap::BFS(vidType source, weight_type *distances) {
  eidType unexplored_edges = M;
  eidType unvisited_vertices = N;
  Direction dir = Direction::TOP_DOWN;
  add_to_frontier(this_frontier, source);
  eidType edges_frontier = rowptr[source + 1] - rowptr[source];
  vidType vertices_frontier = 1;
  distances[source] = 0;
  weight_type distance = 1;

  do {
    if (dir == Direction::BOTTOM_UP && vertices_frontier < N / BETA) {
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
    for (vidType i = 0; i < N; i++) {
      this_frontier[i] = false;
      if (next_frontier[i] == true) {
        edges_frontier += rowptr[i + 1] - rowptr[i];
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
  for (vidType i = 0; i < N; i++) {
    this_frontier[i] = false;
    visited[i] = false;
  }
}