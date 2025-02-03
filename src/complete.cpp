#include <cassert>
#include <graph.hpp>
#include <omp.h>
#include <vector>
#include <limits>
#include <string>

/*
  NOTE FOR REVIEWERS:
  In out code we assumed that the sum of the number of vertices and the number 
  of edges of the graphs is smaller that 2^32 (approx. 4.2 billion). This 
  condition is met by all the graphs used in the competition. If the code 
  is tested on even larger graphs, please uncomment  the #define USE_64_BIT 
  in the line after this comment. This ensures that  the IDs of the vertices 
  and edges don't overflow, although it comes at a small performance penalty 
  of using wider arrays.

  Thank you for the lively competition!
*/

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
namespace complete {
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
  inline void set_distance(mergedType i, weight_type distance) {
    merged[i+1] = distance;
    merged[i] = merged[i] | VISITED_MASK;
  }

  void compute_distances(weight_type *distances, vidType source) const {
    #pragma omp parallel for simd schedule(static)
    for (vidType i = 0; i < N; i++) {
      distances[i] = merged[rowptr[i] + 1];
    }
    distances[source] = 0;
  }

  inline void add_to_frontier(frontier &frontier, mergedType v) const {
    frontier.push_back(v);
  }

  #pragma omp declare reduction(vec_add : std::vector<mergedType> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void top_down_step(const frontier &this_frontier, frontier &next_frontier,
                     const weight_type &distance) {
    #pragma omp parallel for reduction(vec_add : next_frontier) schedule(static) if (this_frontier.size() > 50)
    for (const auto &v : this_frontier) {
      mergedType end = v+2+copy_unmarked(v);
      #pragma omp simd
      for (mergedType i = v + 2; i < end; i++) {
        mergedType neighbor = merged[i];
        if (!IS_VISITED(neighbor)) {
          if (copy_unmarked(neighbor) != 1) {
            add_to_frontier(next_frontier, neighbor);
          }
          set_distance(neighbor, distance);
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) override {
    frontier this_frontier;
    mergedType start = rowptr[source];
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, start);
    set_distance(start, 0);
    weight_type distance = 1;
    while (!this_frontier.empty()) {
      frontier next_frontier;
      next_frontier.reserve(this_frontier.size());
      top_down_step(this_frontier, next_frontier, distance);
      distance++;
      this_frontier = std::move(next_frontier);
    }
    compute_distances(distances, source);
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
  inline void set_distance(vidType i, weight_type distance, weight_type *distances) {
    distances[i] = distance;
    visited[i] = true;
  }

  inline void add_to_frontier(frontier &frontier, vidType v, vidType &edges_frontier) {
    frontier.push_back(v);
    edges_frontier += rowptr[v + 1] - rowptr[v];
  }

  #pragma omp declare reduction(vec_add : std::vector<vidType>, std::vector<std::pair<vidType, bool>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void bottom_up_step(frontier this_frontier, frontier &next_frontier,
                      weight_type distance, weight_type *distances, vidType &edges_frontier) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(static)
    for (vidType i = 0; i < N; i++) {
      if (!visited[i]) {
        for (vidType j = rowptr[i]; j < rowptr[i + 1]; j++) {
            if (visited[col[j]] && distances[col[j]] == distance - 1) {
            // If neighbor is in frontier, add this vertex to next frontier
            if (rowptr[i + 1] - rowptr[i] > 1) {
                add_to_frontier(next_frontier, i, edges_frontier);
            }
            set_distance(i, distance, distances);
            break;
            }
        }
      }
    }
  }

  void top_down_step(frontier this_frontier, frontier &next_frontier,
                     weight_type &distance, weight_type *distances, vidType &edges_frontier) {
    #pragma omp parallel for reduction(vec_add : next_frontier)                    \
    reduction(+ : edges_frontier) schedule(static) if (edges_frontier_old > 150)
    for (const auto &v : this_frontier) {
      for (vidType i = rowptr[v]; i < rowptr[v+1]; i++) {
        vidType neighbor = col[i];
        if (!visited[neighbor]) {
          if (rowptr[neighbor + 1] - rowptr[neighbor] > 1) {
            add_to_frontier(next_frontier, neighbor, edges_frontier);
          }
          set_distance(neighbor, distance, distances);
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) {
    frontier this_frontier;
    dir = Direction::TOP_DOWN;
    vidType edges_frontier = 0;
    add_to_frontier(this_frontier, source, edges_frontier);
    set_distance(source, 0, distances);
    weight_type distance = 1;
    while (!this_frontier.empty()) {
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
        top_down_step(this_frontier, next_frontier, distance, distances, edges_frontier);
      } else {
        bottom_up_step(this_frontier, next_frontier, distance, distances, edges_frontier);
      }
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

  void bottom_up_step(const bool *this_frontier, bool *next_frontier) {
    #pragma omp parallel for schedule(static)
    for (vidType i = 0; i < N; i++) {
      if (!is_visited(i)) {
        for (mergedType j = rowptr[i]; j < rowptr[i+1]; j++) {
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

  void top_down_step(const bool *this_frontier, bool *next_frontier) {
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < N; v++) {
      if (this_frontier[v] == true) {
        mergedType end = rowptr[v + 1];
        #pragma omp simd
        for (mergedType i = rowptr[v]; i < end; i++) {
          vidType neighbor = col[i];
          if (!is_visited(neighbor)) {
            add_to_frontier(next_frontier, neighbor);
          }
        }
      }
    }
  }

  void BFS(vidType source, weight_type *distances) override {
    dir = Direction::TOP_DOWN;
    add_to_frontier(this_frontier, source);
    mergedType edges_frontier = rowptr[source + 1] - rowptr[source];
    vidType vertices_frontier = 1;
    distances[source] = 0;
    weight_type distance = 1;

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
        top_down_step(this_frontier, next_frontier);
      } else {
        bottom_up_step(this_frontier, next_frontier);
      }
        #pragma omp parallel for reduction(+ : edges_frontier, vertices_frontier) schedule(static)
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

inline vidType get_degree(eidType *rowptr, vidType i) {
  return rowptr[i + 1] - rowptr[i];
}

void merged_csr(eidType *rowptr, vidType *col, mergedType *merged, uint64_t N,
                uint64_t M) {
  vidType merged_index = 0;
  // Add degree of each vertex to the start of its neighbor list
  for (vidType i = 0; i < N; i++) {
    vidType start = rowptr[i];
    merged[merged_index++] = get_degree(rowptr, i);
    merged[merged_index++] = std::numeric_limits<mergedType>::max();
    for (vidType j = start; j < rowptr[i + 1]; j++, merged_index++) {
      merged[merged_index] = rowptr[col[j]] + 2*col[j];
    }
  }
  // Fix rowptr indices by adding offset caused by adding the degree to the
  // start of each neighbor list
  for (vidType i = 0; i <= N; i++) {
    rowptr[i] = rowptr[i] + 2*i;
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M, std::string algorithm) {
  bool is_large = algorithm == "large";
  bool is_small = algorithm == "small";
  bool is_classic = algorithm == "classic";

  if (((float)M/N < 10 && !is_small && !is_classic) || is_large) {
    mergedType *merged = new mergedType[M + 2*N];
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