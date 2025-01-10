#include <cassert>
#include <chrono>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <profiling.hpp>

namespace bfs_hybrid_bitmap {

#define ALPHA 14
#define BETA 24
#define VISITED_BIT 31
#define MARKED 1 << VISITED_BIT

enum class Direction { TOP_DOWN, BOTTOM_UP };

class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  uint32_t unexplored_edges, edges_frontier, vertices_frontier;
  Direction dir;
  bool *this_frontier, *next_frontier, *visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, bool *this_frontier, bool *next_frontier, bool *visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), this_frontier(this_frontier), next_frontier(next_frontier), visited(visited), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}

  inline bool is_visited(vidType i) { return visited[i]; }

  inline bool is_unconnected(vidType i) { return rowptr[i] == rowptr[i+1]; }

  void print_frontier(bool* frontier) {
    std::cout << "Frontier: ";
    for (vidType i = 0; i < N; i++) {
      if (frontier[i]) {
        std::cout << i << " ";
      }
    }
    std::cout << std::endl;
  }

  inline void add_to_frontier(bool *frontier, weight_type *distances, vidType v, weight_type distance) {
    frontier[v] = true;
    visited[v] = true;
    vertices_frontier += 1;
    // edges_frontier += rowptr[v+1] - rowptr[v];
    // distances[v] = distance;
  }
  
  void bottom_up_step(bool* this_frontier, bool* next_frontier, weight_type distance, weight_type *distances) {
    // std::cout << "Bottom up step\n";
    #pragma omp parallel for reduction(+:edges_frontier, vertices_frontier) schedule(static, 128)
    for (vidType i = 0; i < N; i++) {
      if (is_visited(i)) {
        continue;
      }
      for (vidType j = rowptr[i]; j < rowptr[i+1]; j++) {
        vidType neighbor = col[j];
        if (this_frontier[neighbor] == true) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, distances, i, distance);
          break;
        }
      }
    }

    #pragma omp parallel for schedule(static, 128)
    for (int i=0 ; i<N; i++){
      if (next_frontier[i] == true) {
        edges_frontier += rowptr[i+1] - rowptr[i];
        vertices_frontier += 1;
        distances[i] = distance;
      }
    }
  }

  void top_down_step(bool* this_frontier, bool* next_frontier, weight_type distance, weight_type *distances) {
    // std::cout << "Top down step\n";
    #pragma omp parallel for reduction(+:edges_frontier, vertices_frontier) schedule(static, 128)
    for (int v = 0; v < N; v++) {
      if (this_frontier[v] == true) {
        for (vidType i = rowptr[v]; i < rowptr[v + 1]; i++) {
          vidType neighbor = col[i];
          if (!is_visited(neighbor)) {
            add_to_frontier(next_frontier, distances, neighbor, distance);
          }
        }
      }
    }

    #pragma omp parallel for schedule(static, 128)
    for (int i=0 ; i<N; i++){
      if (next_frontier[i] == true) {
        edges_frontier += rowptr[i+1] - rowptr[i];
        vertices_frontier += 1;
        distances[i] = distance;
      }
    }
  }

  void BFS(vidType source, weight_type *distances) {
    auto t1 = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("BFS");
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, distances, source, 0);
    edges_frontier += rowptr[source+1] - rowptr[source];
    vertices_frontier += 1;
    distances[source] = 0;
    weight_type distance = 1;
    do {
      if (dir == Direction::BOTTOM_UP && vertices_frontier < N / BETA) {
        dir = Direction::TOP_DOWN;
      } else if (dir == Direction::TOP_DOWN && edges_frontier > unexplored_edges / ALPHA) {
        dir = Direction::BOTTOM_UP;
      }
      unexplored_edges -= edges_frontier;
      edges_frontier = 0;
      vertices_frontier = 0;
      if (dir == Direction::TOP_DOWN) {
        top_down_step(this_frontier, next_frontier, distance, distances);
      } else {
        bottom_up_step(this_frontier, next_frontier, distance, distances);
      }
      if (vertices_frontier == 0) {
        break;
      }
      std::swap(this_frontier, next_frontier);
      #pragma omp parallel for
      for (int i=0 ; i<N; i++){
        next_frontier[i] = false;
      }
      distance++;
    } while(true);
    LIKWID_MARKER_STOP("BFS");
    LIKWID_MARKER_CLOSE;
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "BFS: " << ms_double.count() << "ms\n";
  }
};

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  auto t1 = std::chrono::high_resolution_clock::now();
  bool *this_frontier = new bool[N];
  bool *next_frontier = new bool[N];
  bool *visited = new bool[N];
  std::fill(this_frontier, this_frontier + N, false);
  std::fill(next_frontier, next_frontier + N, false);
  std::fill(visited, visited + N, false);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Preprocessing: " << ms_double.count() << "ms\n";
  return new Graph(rowptr, col, this_frontier, next_frontier, visited, N, M);
}

}