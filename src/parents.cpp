#include <cassert>
#include <graph.hpp>
#include <omp.h>
#include <vector>
#include <string>
//#define USE_64_BIT

#ifdef USE_64_BIT
using mergedType = uint64_t;
constexpr int VISITED_BIT = 63;
constexpr long int VISITED_MASK = 1L << VISITED_BIT;
#else
using mergedType = uint32_t;
constexpr int VISITED_BIT = 31;
constexpr int VISITED_MASK = 1 << VISITED_BIT;
#endif

constexpr int ALPHA = 4;
constexpr int BETA = 24;

enum class Direction { TOP_DOWN, BOTTOM_UP };
using frontier = std::vector<mergedType>;

#define IS_VISITED(i) (((merged[i]) & VISITED_MASK) != 0)
namespace parents {
namespace large_graph {
class Graph final : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  mergedType *merged;
  Direction dir;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, mergedType *merged, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), merged(merged), N(N), M(M) {}
  ~Graph() = default;
  inline mergedType copy_unmarked(mergedType i) const {
    return merged[i] & ~VISITED_MASK;
  }
  inline void set_parent(mergedType i, weight_type parent) {
    merged[i+1] = parent;
    merged[i] = merged[i] | VISITED_MASK;
  }

  void compute_parents(weight_type *parents, vidType source) const {
    #pragma omp parallel for simd schedule(static)
    for (vidType i = 0; i < N; i++) {
      parents[i] = merged[rowptr[i] + 1];
    }
  }

  inline void add_to_frontier(frontier &frontier, mergedType v) const {
    frontier.push_back(v);
  }

  #pragma omp declare reduction(vec_add : std::vector<mergedType> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void top_down_step(const frontier &this_frontier, frontier &next_frontier) {
    #pragma omp parallel for reduction(vec_add : next_frontier) schedule(static) if (this_frontier.size() > 50)
    for (const auto &v : this_frontier) {
      vidType end = v + merged[v+2] + 3;
      #pragma omp simd
      for (mergedType i = v + 3; i < end; i++) {
        mergedType neighbor = merged[i];
        if (!IS_VISITED(neighbor)) {
          add_to_frontier(next_frontier, neighbor);
          set_parent(neighbor, copy_unmarked(v));
        }
      }
    }
  }

  void BFS(vidType source, weight_type *parents) override {
    frontier this_frontier;
    mergedType start = rowptr[source];
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, start);
    set_parent(start, source);
    while (!this_frontier.empty()) {
      frontier next_frontier;
      next_frontier.reserve(this_frontier.size());
      top_down_step(this_frontier, next_frontier);
      this_frontier = std::move(next_frontier);
    }
    compute_parents(parents, source);
  }
};
} // namespace large_graph

namespace classic {
class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  uint32_t unexplored_edges;
  Direction dir;
  bool* visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;
  uint32_t edges_frontier_old;

public:
  Graph(eidType *rowptr, vidType *col, bool* visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), visited(visited), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}
  inline void set_parent(vidType i, weight_type parent, weight_type *parents) {
    parents[i] = parent;
    visited[i] = true;
  }

  inline void add_to_frontier(frontier &frontier, vidType v) {
    frontier.push_back(v);
  }

  #pragma omp declare reduction(vec_add : std::vector<vidType>, std::vector<std::pair<vidType, bool>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void top_down_step(frontier this_frontier, frontier &next_frontier,
                     weight_type &distance, weight_type *parents) {
    #pragma omp parallel for reduction(vec_add : next_frontier) schedule(static) if (edges_frontier_old > 150)
    for (const auto &v : this_frontier) {
      for (vidType i = rowptr[v]; i < rowptr[v+1]; i++) {
        vidType neighbor = col[i];
        if (!visited[neighbor]) {
          if (rowptr[neighbor + 1] - rowptr[neighbor] > 1) {
            add_to_frontier(next_frontier, neighbor);
          }
          set_parent(neighbor, v, parents);
        }
      }
    }
  }

  void BFS(vidType source, weight_type *parents) {
    frontier this_frontier;
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, source);
    set_parent(source, source, parents);
    weight_type distance = 1;
    while (!this_frontier.empty()) {
      // std::cout << "Edges in frontier " << edges_frontier << ", vertices in frontier " << this_frontier.size() << ", unexplored edges " << unexplored_edges << std::endl;
      frontier next_frontier;
      next_frontier.reserve(this_frontier.size());
      top_down_step(this_frontier, next_frontier, distance, parents);
      distance++;
      this_frontier = std::move(next_frontier);
    }
  }
};
} // namespace classic

namespace small_graph {

class Graph : public BaseGraph {
  mergedType *rowptr;
  [[maybe_unused]] vidType *col;
  bool *this_frontier;
  bool *next_frontier;
  bool *visited;
  mergedType unexplored_edges;
  mergedType unvisited_vertices;
  Direction dir;
  std::vector<vidType> to_visit;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(mergedType *rowptr, vidType *col, bool *this_frontier, bool *next_frontier,
        bool *visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), this_frontier(this_frontier),
        next_frontier(next_frontier), visited(visited), unexplored_edges(M), unvisited_vertices(N), N(N), M(M) {}
  ~Graph() = default;

  inline bool is_visited(vidType i) const { return visited[i]; }

  inline void add_to_frontier(bool *frontier, vidType v) {
    frontier[v] = true; //__builtin_nontemporal_store(true, &frontier[v]);
    visited[v] = true; //__builtin_nontemporal_store(true, &visited[v]);
  }

  void bottom_up_step(const bool *this_frontier, bool *next_frontier, weight_type *parents) {
    #pragma omp parallel for schedule(static)
    for (vidType i = 0; i < N; i++) {
      if (!is_visited(i)) {
        for (mergedType j = rowptr[i]; j < rowptr[i+1]; j++) {
          vidType neighbor = col[j];
          if (this_frontier[neighbor] == true) {
            // If neighbor is in frontier, add this vertex to next frontier
            add_to_frontier(next_frontier, i);
            parents[i] = neighbor;
            break;
          }
        }
      }
    }
  }

  void top_down_step(const bool *this_frontier, bool *next_frontier, weight_type *parents) {
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < N; v++) {
      if (this_frontier[v] == true) {
        mergedType end = rowptr[v + 1];
        #pragma omp simd
        for (mergedType i = rowptr[v]; i < end; i++) {
          vidType neighbor = col[i];
          if (!is_visited(neighbor)) {
            add_to_frontier(next_frontier, neighbor);
            parents[neighbor] = v;
          }
        }
      }
    }
  }

  void BFS(vidType source, weight_type *parents) override {
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, source);
    mergedType edges_frontier = rowptr[source + 1] - rowptr[source];
    vidType vertices_frontier = 1;
    parents[source] = source;

    do {
      if (dir == Direction::BOTTOM_UP && vertices_frontier < N / BETA) {
        dir = Direction::TOP_DOWN;
      } else if (dir == Direction::TOP_DOWN && edges_frontier > unexplored_edges / ALPHA) {
        dir = Direction::BOTTOM_UP;
      }
      unexplored_edges -= edges_frontier;
      unvisited_vertices -= vertices_frontier;
      edges_frontier = 0;
      vertices_frontier = 0;
      if (dir == Direction::TOP_DOWN) {
        top_down_step(this_frontier, next_frontier, parents);
      } else {
        bottom_up_step(this_frontier, next_frontier, parents);
      }
        #pragma omp parallel for reduction(+ : edges_frontier, vertices_frontier) schedule(static)
        for (vidType i = 0; i < N; i++) {
          this_frontier[i] = false;
          if (next_frontier[i] == true) {
            edges_frontier += rowptr[i + 1] - rowptr[i];
            vertices_frontier += 1;
          }
        }
      if (vertices_frontier == 0) {
        break;
      }
      std::swap(this_frontier, next_frontier);
    } while (true);
  }
};
} // namespace small_graph

inline vidType get_degree(eidType *rowptr, vidType i) {
  return rowptr[i + 1] - rowptr[i];
}

void merged_csr(eidType *rowptr, vidType *col, mergedType *merged, uint64_t N,
                uint64_t M) {
  vidType merged_index = 0;
  // Add degree of each vertex to the start of its neighbor list
  for (vidType i = 0; i < N; i++) {
    vidType start = rowptr[i];
    merged[merged_index++] = i;
    merged[merged_index++] = -1;
    merged[merged_index++] = get_degree(rowptr, i);
    for (vidType j = start; j < rowptr[i + 1]; j++, merged_index++) {
      merged[merged_index] = rowptr[col[j]] + 3*col[j];
    }
  }
  // Fix rowptr indices by adding offset caused by adding the degree to the
  // start of each neighbor list
  for (vidType i = 0; i <= N; i++) {
    rowptr[i] = rowptr[i] + 3*i;
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M, std::string algorithm) {
  bool is_large = algorithm == "large";
  bool is_small = algorithm == "small";
  bool is_classic = algorithm == "classic";

  if (((float)M/N < 10 && !is_small && !is_classic) || is_large) {
    mergedType *merged = new mergedType[M + 3*N];
    merged_csr(rowptr, col, merged, N, M);
    return new large_graph::Graph(rowptr, col, merged, N, M);
  } else if (!is_classic && !is_large) {
    bool *this_frontier = new bool[N];
    mergedType *newrowptr = new mergedType[N+1];
    bool *next_frontier = new bool[N];
    bool *visited = new bool[N];
    #pragma omp parallel for schedule(static)
    for (mergedType i = 0; i < N; i++) {
      newrowptr[i] = rowptr[i];
      this_frontier[i] = false;
      next_frontier[i] = false;
      visited[i] = false;
    }
    newrowptr[N] = rowptr[N];
    return new small_graph::Graph(newrowptr, col, this_frontier, next_frontier, visited, N, M);
  } else {
    bool *visited = new bool[N];
    #pragma omp parallel for schedule(static)
    for (mergedType i = 0; i < N; i++) {
      visited[i] = false;
    }
    return new classic::Graph(rowptr, col, visited, N, M);
  }
}

}