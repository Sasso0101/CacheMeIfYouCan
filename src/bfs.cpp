#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <vector>
#include <profiling.hpp>

namespace bfs {

#define VISITED_BIT 30
#define MARKED_BIT 31
#define MARKED 1 << MARKED_BIT
#define VISITED 0b11 << VISITED_BIT
#define VISITED_MASK 0b01 << VISITED_BIT
#define ALPHA 14
#define BETA 24

enum class Direction { TOP_DOWN, BOTTOM_UP };
typedef std::vector<vidType> frontier;

class Graph : public BaseGraph {
  eidType *rowptr;
  [[maybe_unused]] vidType *col;
  vidType *merged;
  uint32_t unexplored_edges;
  Direction dir;
  bool *visited;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;
  uint32_t edges_frontier;

public:
  Graph(eidType *rowptr, vidType *col, vidType *merged, bool *visited, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), merged(merged), visited(visited), N(N), M(M) {
    unexplored_edges = M;
  }
  ~Graph() {}

  inline vidType copy_unmarked(vidType i) { return merged[i] & ~(VISITED); }

  inline bool is_marked(vidType i) { return ((merged[i]) & MARKED) != 0; }
  // inline bool is_visited(vidType i) { return ((merged[i] & VISITED_MASK) >> VISITED_BIT) != 0; }
  inline bool is_visited(vidType i) { return visited[i]; }

  inline void set_distance(vidType i, weight_type distance) {
    merged[i] = distance | VISITED;
    visited[i] = true;
  }

  inline bool is_unconnected(vidType i) { return (merged[i] & ~(VISITED)) == 0; }

  inline bool is_leaf(vidType i) { return (rowptr[i] == rowptr[i + 1] - 1); }

  void compute_distances(weight_type *distances, vidType source) {
    #pragma omp parallel for schedule(auto)
    for (uint64_t i = 0; i < N; i++) {
      if (is_visited(rowptr[i])) {
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

  bool switch_to_bottom_up(frontier &frontier) {
    bool to_switch = edges_frontier > unexplored_edges / ALPHA;
    unexplored_edges -= edges_frontier;
    // std::cout << "Edges in frontier: " << edges_frontier << ", unexplored edges: " << unexplored_edges << ", switch: " << to_switch << std::endl;
    edges_frontier = 0;
    return to_switch;
  }

  bool switch_to_top_down(frontier &frontier) {
    bool to_switch = frontier.size() < N / BETA;
    unexplored_edges -= edges_frontier;
    // std::cout << "Frontier size: " << frontier->size() << ", total elements: " << N << ", switch: " << to_switch << std::endl;
    edges_frontier = 0;
    return to_switch;
  }

  void bottom_up_step(frontier this_frontier, frontier &next_frontier, weight_type distance) {
    // std::cout << "Bottom up step\n";
    #pragma omp parallel for reduction(vec_add: next_frontier) reduction(+: edges_frontier) schedule(auto)
    for (vidType i = 0; i < N; i++) {
      vidType start = rowptr[i];
      if (is_visited(start) || is_unconnected(start)) {
        continue;
      }
      for (vidType j = start + 1; j < rowptr[i+1]; j++) {
        if (is_visited(merged[j]) && copy_unmarked(merged[j]) == distance - 1) {
          // If neighbor is in frontier, add this vertex to next frontier
          add_to_frontier(next_frontier, start);
          set_distance(start, distance);
          break;
        }
      }
    }
  }

  #pragma omp declare reduction(vec_add: std::vector<vidType>: \
    omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  void top_down_step(frontier this_frontier, frontier &next_frontier, weight_type distance) {
    // std::cout << "Top down step\n";
    #pragma omp parallel for reduction(vec_add: next_frontier) reduction(+: edges_frontier) schedule(auto)
    for (const auto &v : this_frontier) {
      for (vidType i = v + 1; !is_marked(i) && (i < N + M); ++i) {
        vidType neighbor = merged[i];
        if (!is_visited(neighbor)) {
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
      frontier next_frontier;
      next_frontier.reserve(this_frontier.size());
      if (dir == Direction::TOP_DOWN) {
        if (switch_to_bottom_up(this_frontier)) {
          dir = Direction::BOTTOM_UP;
          bottom_up_step(this_frontier, next_frontier, distance);
        } else {
          top_down_step(this_frontier, next_frontier, distance);
        }
      } else {
        if (dir == Direction::BOTTOM_UP) {
          if (switch_to_top_down(this_frontier)) {
          dir = Direction::TOP_DOWN;
          top_down_step(this_frontier, next_frontier, distance);
          } else {
            bottom_up_step(this_frontier, next_frontier, distance);
          }
        }
      }
      distance++;
      this_frontier = std::move(next_frontier);
    }
    compute_distances(distances, source);
  }
};

inline vidType get_degree(eidType *rowptr, vidType i) {
  return rowptr[i + 1] - rowptr[i];
}

void merged_csr(eidType *rowptr, vidType *col, vidType *merged, uint64_t N, uint64_t M) {
  vidType merged_index = 0;
  // Add degree of each vertex to the start of its neighbor list
  for (vidType i = 0; i < N; i++) {
    vidType start = rowptr[i];
    merged[merged_index++] = get_degree(rowptr, i) | MARKED;
    for (vidType j = start; j < rowptr[i + 1]; j++, merged_index++) {
      merged[merged_index] = rowptr[col[j]] + col[j];
    }
  }
  // Fix rowptr indices by adding offset caused by adding the degree to the start of
  // each neighbor list
  for (vidType i = 0; i <= N; i++) {
    rowptr[i] = rowptr[i] + i;
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  vidType *merged = new vidType[M+N];
  merged_csr(rowptr, col, merged, N, M);
  bool *visited = new bool[N+M];
  std::fill(visited, visited + N + M, false);
  return new Graph(rowptr, col, merged, visited, N, M);
}

} // namespace bfs
