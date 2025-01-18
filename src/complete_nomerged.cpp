#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <profiling.hpp>
#include <vector>

constexpr int ALPHA = 3;
constexpr int BETA = 12;

enum class Direction { TOP_DOWN, BOTTOM_UP };
typedef std::vector<vidType> frontier;

namespace large_graph {
class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  uint32_t unexplored_edges;
  Direction dir;
  bool* visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;
  uint32_t edges_frontier, edges_frontier_old;

public:
  Graph(eidType *rowptr, vidType *col, bool* visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), visited(visited), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}
  inline void set_distance(vidType i, weight_type distance, weight_type *distances) {
    distances[i] = distance;
    visited[i] = true;
  }

  void print_frontier(frontier &frontier) {
    std::cout << "Frontier: ";
    for (const auto &v : frontier) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  inline void add_to_frontier(frontier &frontier, vidType v) {
    frontier.push_back(v);
    edges_frontier += rowptr[v + 1] - rowptr[v];
  }

  #pragma omp declare reduction(vec_add : std::vector<vidType>, std::vector<std::pair<vidType, bool>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void bottom_up_step(frontier this_frontier, frontier &next_frontier,
                      weight_type distance, weight_type *distances) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(auto)
    for (vidType i = 0; i < N; i++) {
      if (visited[i]) {
        continue;
      }
      for (vidType j = rowptr[i]; j < rowptr[i + 1]; j++) {
        if (visited[col[j]] && distances[col[j]] == distance - 1) {
          // If neighbor is in frontier, add this vertex to next frontier
          if (rowptr[i + 1] - rowptr[i] > 1) {
            add_to_frontier(next_frontier, i);
          }
          set_distance(i, distance, distances);
          break;
        }
      }
    }
  }

  void top_down_step(frontier this_frontier, frontier &next_frontier,
                     weight_type &distance, weight_type *distances) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(auto) if (edges_frontier_old > 150)
    for (const auto &v : this_frontier) {
      for (vidType i = rowptr[v]; i < rowptr[v+1]; i++) {
        vidType neighbor = col[i];
        if (!visited[neighbor]) {
          if (rowptr[neighbor + 1] - rowptr[neighbor] > 1) {
            add_to_frontier(next_frontier, neighbor);
          }
          set_distance(neighbor, distance, distances);
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) {
    frontier this_frontier;
    vidType start = rowptr[source];
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, start);
    set_distance(start, 0, distances);
    weight_type distance = 1;
    while (!this_frontier.empty()) {
      // print_frontier(this_frontier);
      // std::cout << "Edges in frontier " << edges_frontier << ", vertices in frontier " << this_frontier.size() << ", unexplored edges " << unexplored_edges << std::endl;
      frontier next_frontier;
      next_frontier.reserve(this_frontier.size());
      if (dir == Direction::BOTTOM_UP && this_frontier.size() < N / BETA) {
        dir = Direction::TOP_DOWN;
      } else if (dir == Direction::TOP_DOWN &&
                 edges_frontier > unexplored_edges / ALPHA) {
        dir = Direction::BOTTOM_UP;
      }
      unexplored_edges -= edges_frontier;
      edges_frontier_old = edges_frontier;
      edges_frontier = 0;
      if (dir == Direction::TOP_DOWN) {
        top_down_step(this_frontier, next_frontier, distance, distances);
      } else {
        bottom_up_step(this_frontier, next_frontier, distance, distances);
      }
      distance++;
      this_frontier = std::move(next_frontier);
    }
  }
};
} // namespace large_graph

namespace small_graph {
class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  uint32_t unexplored_edges, edges_frontier, vertices_frontier,
      unvisited_vertices;
  Direction dir;
  bool *this_frontier, *next_frontier, *visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, bool *this_frontier, bool *next_frontier,
        bool *visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), this_frontier(this_frontier),
        next_frontier(next_frontier), visited(visited), N(N), M(M) {
    unexplored_edges = M;
    unvisited_vertices = N;
  }
  ~Graph() {}

  inline bool is_visited(vidType i) { return visited[i]; }

  inline bool is_unconnected(vidType i) { return rowptr[i] == rowptr[i + 1]; }

  void print_frontier(bool *frontier) {
    std::cout << "Frontier: ";
    for (vidType i = 0; i < N; i++) {
      if (frontier[i]) {
        std::cout << i << " ";
      }
    }
    std::cout << std::endl;
  }

  inline void add_to_frontier(bool *frontier, weight_type *distances, vidType v,
                              weight_type distance) {
    frontier[v] = true;
    visited[v] = true;
  }

  #pragma omp declare reduction(vec_add : std::vector<vidType> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void bottom_up_step(bool *this_frontier, bool *next_frontier,
                      weight_type distance, weight_type *distances) {
    // std::cout << "Bottom up step\n";
    #pragma omp parallel for schedule(auto)
    for (vidType i = 0; i < N; i++) {
      if (is_visited(i)) {
        continue;
      }
      for (vidType j = rowptr[i]; j < rowptr[i + 1]; j++) {
        vidType neighbor = col[j];
        if (this_frontier[neighbor] == true) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, distances, i, distance);
          break;
        }
      }
    }
  }

  void top_down_step(bool *this_frontier, bool *next_frontier,
                     weight_type distance, weight_type *distances) {
    // std::cout << "Top down step\n";
    #pragma omp parallel for schedule(auto)
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
  }

  void BFS(vidType source, weight_type *distances) {
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, distances, source, 0);
    edges_frontier = rowptr[source + 1] - rowptr[source];
    vertices_frontier = 1;
    distances[source] = 0;
    weight_type distance = 1;

    do {
      if (dir == Direction::BOTTOM_UP && vertices_frontier < N / BETA) {
        dir = Direction::TOP_DOWN;
      } else if (dir == Direction::TOP_DOWN &&
                 edges_frontier > unexplored_edges / ALPHA) {
        dir = Direction::BOTTOM_UP;
      }
      if (dir == Direction::TOP_DOWN) {
        top_down_step(this_frontier, next_frontier, distance, distances);
      } else {
        bottom_up_step(this_frontier, next_frontier, distance, distances);
      }
      unexplored_edges -= edges_frontier;
      unvisited_vertices -= vertices_frontier;
      edges_frontier = 0;
      vertices_frontier = 0;

      #pragma omp parallel for reduction(+ : edges_frontier, vertices_frontier) schedule(auto)
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
  }
};
} // namespace small_graph

namespace complete {
  inline vidType get_degree(eidType *rowptr, vidType i) {
    return rowptr[i + 1] - rowptr[i];
  }

  BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                              uint64_t M) {
    if ((float)M/N < 10) {
      bool* visited = new bool[N];
      for (vidType i = 0; i < N; i++) {
        visited[i] = false;
      }
      return new large_graph::Graph(rowptr, col, visited, N, M);
    } else {
      bool *this_frontier = new bool[N];
      bool *next_frontier = new bool[N];
      bool *visited = new bool[N];
      std::fill(this_frontier, this_frontier + N, false);
      std::fill(next_frontier, next_frontier + N, false);
      std::fill(visited, visited + N, false);
      return new small_graph::Graph(rowptr, col, this_frontier, next_frontier,
                                    visited, N, M);
    }
  }
}
