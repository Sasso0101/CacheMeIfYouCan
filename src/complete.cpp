#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <profiling.hpp>
#include <vector>

constexpr int ALPHA = 3;
constexpr int BETA = 12;
constexpr int MARKED_BIT = 31;
constexpr int MARKED_MASK = 1 << MARKED_BIT;
constexpr int VISITED_BIT = 30;
constexpr int VISITED_MASK = 0b01 << VISITED_BIT;
constexpr int MARKED_VISITED_MASK = 0b11 << VISITED_BIT;

enum class Direction { TOP_DOWN, BOTTOM_UP };
typedef std::vector<vidType> frontier;

#define IS_VISITED(i) (((merged[i]) & VISITED_MASK) != 0)
#define IS_MARKED(i) (((merged[i]) & MARKED_MASK) != 0)

namespace large_graph {
class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  vidType *merged;
  uint32_t unexplored_edges;
  Direction dir;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;
  uint32_t edges_frontier, edges_frontier_old;

public:
  Graph(eidType *rowptr, vidType *col, vidType *merged, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), merged(merged), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}
  inline vidType copy_unmarked(vidType i) {
    return merged[i] & ~(MARKED_VISITED_MASK);
  }
  inline void set_distance(vidType i, weight_type distance) {
    merged[i] = distance | MARKED_VISITED_MASK;
  }

  inline bool is_unconnected(vidType i) {
    return (merged[i] & ~(MARKED_VISITED_MASK)) == 0;
  }

  inline bool is_leaf(vidType i) { return (rowptr[i] == rowptr[i + 1] - 1); }

  void compute_distances(weight_type *distances, vidType source) {
    #pragma omp parallel for schedule(auto)
    for (uint64_t i = 0; i < N; i++) {
      if (IS_VISITED(rowptr[i])) {
        distances[i] = copy_unmarked(rowptr[i]);
      }
    }
    distances[source] = 0;
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
    edges_frontier += copy_unmarked(v);
  }

  #pragma omp declare reduction(vec_add : std::vector<vidType>, std::vector<std::pair<vidType, bool>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void bottom_up_step(frontier this_frontier, frontier &next_frontier,
                      weight_type distance) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(auto)
    for (vidType i = 0; i < N; i++) {
      vidType start = rowptr[i];
      if (IS_VISITED(start)) {
        continue;
      }
      for (vidType j = start + 1; j < rowptr[i + 1]; j++) {
        if (IS_VISITED(merged[j]) && copy_unmarked(merged[j]) == distance - 1) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, start);
          set_distance(start, distance);
          break;
        }
      }
    }
  }

  void top_down_step(frontier this_frontier, frontier &next_frontier,
                     weight_type &distance) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(auto) if (edges_frontier_old > 150)
    for (const auto &v : this_frontier) {
      for (vidType i = v + 1; !IS_MARKED(i); i++) {
        vidType neighbor = merged[i];
        if (!IS_VISITED(neighbor)) {
          add_to_frontier(next_frontier, neighbor);
          set_distance(neighbor, distance);
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) {
    frontier this_frontier;
    vidType start = rowptr[source];
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, start);
    set_distance(start, 0);
    weight_type distance = 1;
    while (!this_frontier.empty()) {
      // print_frontier(this_frontier);
      // std::cout << "Edges in frontier " << edges_frontier << ", vertices in
      // frontier " << this_frontier.size();
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
        top_down_step(this_frontier, next_frontier, distance);
      } else {
        bottom_up_step(this_frontier, next_frontier, distance);
      }
      distance++;
      this_frontier = std::move(next_frontier);
    }
    compute_distances(distances, source);
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

  void merged_csr(eidType *rowptr, vidType *col, vidType *merged, uint64_t N,
                  uint64_t M) {
    vidType merged_index = 0;
    // Add degree of each vertex to the start of its neighbor list
    for (vidType i = 0; i < N; i++) {
      vidType start = rowptr[i];
      merged[merged_index++] = get_degree(rowptr, i) | MARKED_MASK;
      for (vidType j = start; j < rowptr[i + 1]; j++, merged_index++) {
        merged[merged_index] = rowptr[col[j]] + col[j];
      }
    }
    // Fix rowptr indices by adding offset caused by adding the degree to the
    // start of each neighbor list
    for (vidType i = 0; i <= N; i++) {
      rowptr[i] = rowptr[i] + i;
    }
    merged[M + N] = MARKED_MASK;
  }

  BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                              uint64_t M) {
    if ((float)M/N < 10) {
      vidType *merged = new vidType[M + N + 1];
      merged_csr(rowptr, col, merged, N, M);
      return new large_graph::Graph(rowptr, col, merged, N, M);
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
