#include <chrono>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <algorithm>
// #include <omp.h>
#include <unordered_map>
#include <vector>
#include <profiling.hpp>

namespace bfs_hybrid_bitmap {

#define ALPHA 14
#define BETA 24

typedef std::unordered_map<vidType, vidType> degrees_map;
enum class Direction { TOP_DOWN, BOTTOM_UP };
typedef std::vector<vidType> frontier;

class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  uint32_t unexplored_edges;
  Direction dir;
  bool *visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;
  uint32_t edges_frontier;

public:
  Graph(eidType *rowptr, vidType *col, bool *visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), visited(visited), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}

  inline bool is_visited(vidType i) { return visited[i]; }

  inline bool is_unconnected(vidType i) { return rowptr[i] == rowptr[i+1]; }

  void print_frontier(frontier &frontier) {
    std::cout << "Frontier: ";
    for (const auto &v : frontier) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  inline void add_to_frontier(frontier *frontier, vidType v) {
    frontier->push_back(v);
    edges_frontier += rowptr[v + 1] - rowptr[v];
  }

  bool switch_to_bottom_up(frontier *frontier) {
    bool to_switch = edges_frontier > unexplored_edges / ALPHA;
    unexplored_edges -= edges_frontier;
    // std::cout << "Edges in frontier: " << edges_frontier << ", unexplored edges: " << unexplored_edges << ", switch: " << to_switch << std::endl;
    edges_frontier = 0;
    return to_switch;
  }

  bool switch_to_top_down(frontier *frontier) {
    bool to_switch = frontier->size() < N / BETA;
    unexplored_edges -= edges_frontier;
    // std::cout << "Frontier size: " << frontier->size() << ", total elements: " << N << ", switch: " << to_switch << std::endl;
    edges_frontier = 0;
    return to_switch;
  }

  void bottom_up_step(frontier this_frontier, frontier *next_frontier, weight_type distance, weight_type *distances) {
    // std::cout << "Bottom up step\n";
    for (vidType i = 0; i < N; i++) {
      vidType start = rowptr[i];
      if (is_visited(i) || is_unconnected(i)) {
        continue;
      }
      for (vidType j = start; j < rowptr[i+1]; j++) {
        vidType neighbor = col[j];
        if (is_visited(neighbor) && distances[neighbor] == distance - 1) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, i);
          distances[i] = distance;
          visited[i] = true;
          break;
        }
      }
    }
  }

  void top_down_step(frontier this_frontier, frontier *next_frontier, weight_type distance, weight_type *distances) {
    // std::cout << "Top down step\n";
    for (const auto &v : this_frontier) {
      for (vidType i = rowptr[v]; i < rowptr[v + 1]; i++) {
        vidType neighbor = col[i];
        if (!is_visited(neighbor)) {
          add_to_frontier(next_frontier, neighbor);
          distances[neighbor] = distance;
          visited[neighbor] = true;
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) {
    auto t1 = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("BFS");
    frontier this_frontier;
    dir = Direction::TOP_DOWN;
    add_to_frontier(&this_frontier, source);
    distances[source] = 0;
    visited[source] = true;
    weight_type distance = 1;
    while (!this_frontier.empty()) {
      // print_frontier(this_frontier);
      frontier next_frontier;
      if (dir == Direction::TOP_DOWN) {
        if (switch_to_bottom_up(&this_frontier)) {
          dir = Direction::BOTTOM_UP;
          bottom_up_step(this_frontier, &next_frontier, distance, distances);
        } else {
          top_down_step(this_frontier, &next_frontier, distance, distances);
        }
      } else {
        if (dir == Direction::BOTTOM_UP) {
          if (switch_to_top_down(&this_frontier)) {
          dir = Direction::TOP_DOWN;
          top_down_step(this_frontier, &next_frontier, distance, distances);
          } else {
            bottom_up_step(this_frontier, &next_frontier, distance, distances);
          }
        }
      }
      distance++;
      std::swap(this_frontier, next_frontier);
    }
    LIKWID_MARKER_STOP("BFS");
    LIKWID_MARKER_CLOSE;
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "BFS: " << ms_double.count() << "ms\n";
  }
};

inline vidType get_degree(eidType *rowptr, vidType i) {
  return rowptr[i + 1] - rowptr[i];
}


BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  auto t1 = std::chrono::high_resolution_clock::now();
  bool *visited = new bool[N];
  std::fill(visited, visited + N, false);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Preprocessing: " << ms_double.count() << "ms\n";
  return new Graph(rowptr, col, visited, N, M);
}

}