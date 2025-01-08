#include <chrono>
#include <cstdint>
#include <graph.hpp>
#include <iostream>
#include <algorithm>
// #include <omp.h>
#include <unordered_map>
#include <vector>
#include <profiling.hpp>

namespace bfs {

#define MARKED 1 << 31
#define MSB 31
#define ALPHA 14
#define BETA 24

typedef std::unordered_map<vidType, vidType> degrees_map;
enum class Direction { TOP_DOWN, BOTTOM_UP };

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  degrees_map *degrees;
  uint32_t unexplored_edges;
  Direction dir;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, degrees_map *degrees, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), degrees(degrees), N(N), M(M) {
    unexplored_edges = N;
  }
  ~Graph() {}

  inline vidType copy_unmarked(vidType i) { return col[i] & ~(MARKED); }

  inline bool is_marked(vidType i) { return (col[i] >> MSB) != 0; }

  inline void mark(vidType i) { col[i] = col[i] | MARKED; }

  inline void set_distance(vidType i, weight_type distance) {
    col[i] = distance | MARKED;
  }

  inline bool is_unconnected(vidType i) { return (rowptr[i] == rowptr[i + 1]); }

  inline bool is_leaf(vidType i) { return (rowptr[i] == rowptr[i + 1] - 1); }

  void compute_distances(weight_type *distances, vidType source) {
    for (uint64_t i = 0; i < N; i++) {
      // Check if node has only one neighbor
      if (is_unconnected(i)) {
        distances[i] = std::numeric_limits<weight_type>::max();
      } else if (is_leaf(i)) {
        // For leaf nodes, the neighbor is the ID of the vertex
        vidType neighbour = copy_unmarked(rowptr[i]);
        if (!is_leaf(neighbour) && is_marked(rowptr[neighbour])) {
          distances[i] = copy_unmarked(rowptr[neighbour]) + 1;
        }
      } else if (is_marked(rowptr[i])) {
        distances[i] = copy_unmarked(rowptr[i]);
      }
    }
    distances[source] = 0;
  }

  void print_frontier(std::vector<vidType> &frontier) {
    std::cout << "Frontier: ";
    for (const auto &v : frontier) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }

  uint32_t count_edges_frontier(std::vector<vidType> *frontier) {
    uint32_t edges_frontier = 0;
    for (const auto &v : *frontier) {
      edges_frontier += degrees->at(v);
    }
    return edges_frontier;
  }

  bool switch_to_bottom_up(std::vector<vidType> *frontier) {
    uint32_t edges_frontier = count_edges_frontier(frontier);
    bool to_switch = edges_frontier > unexplored_edges / ALPHA;
    unexplored_edges -= edges_frontier;
    std::cout << "Edges in frontier: " << edges_frontier << ", unexplored edges: " << unexplored_edges << ", switch: " << to_switch << std::endl;
    return to_switch;
  }

  bool switch_to_top_down(std::vector<vidType> *frontier) {
    uint32_t edges_frontier = count_edges_frontier(frontier);
    bool to_switch = frontier->size() < N / BETA;
    unexplored_edges -= edges_frontier;
    std::cout << "Frontier size: " << frontier->size() << ", total elements: " << N << ", switch: " << to_switch << std::endl;
    return to_switch;
  }

  void bottom_up_step(std::vector<vidType> this_frontier, std::vector<vidType> *next_frontier, weight_type distance) {
    std::cout << "Bottom up step\n";
    for (vidType i = 0; i < N; i++) {
      if (is_marked(rowptr[i])) {
        continue;
      }
      vidType start = rowptr[i];
      vidType end = rowptr[i + 1];
      for (vidType j = start; j < end; j++) {
        vidType neighbor = copy_unmarked(j);
        if (!is_marked(neighbor)) {
          continue;
        }
        if (std::find(this_frontier.begin(), this_frontier.end(), neighbor) !=
            this_frontier.end()) {
          next_frontier->push_back(rowptr[i]);
          mark(rowptr[i]);
          break;
        }
      }
    }
    for (const auto &v : this_frontier) {
      set_distance(v, distance);
    }
  }

  void top_down_step(std::vector<vidType> this_frontier, std::vector<vidType> *next_frontier, weight_type distance) {
    std::cout << "Top down step\n";
    for (const auto &v : this_frontier) {
      vidType curr_index = v;
      vidType neighbor_start = copy_unmarked(curr_index);
      // Repeat until all neighbors have been visited except last one
      do {
        if (!is_marked(neighbor_start)) {
          mark(neighbor_start);
          next_frontier->push_back(neighbor_start);
        }
        curr_index++;
        neighbor_start = copy_unmarked(curr_index);
      } while (!is_marked(curr_index));
      // Visit last neighbor
      if (!is_marked(neighbor_start)) {
        mark(neighbor_start);
        next_frontier->push_back(neighbor_start);
      }
      set_distance(v, distance);
    }
  }

  void BFS(vidType source, weight_type *distances) {
    auto t1 = std::chrono::high_resolution_clock::now();
    LIKWID_MARKER_START("BFS");
    std::vector<vidType> this_frontier;
    weight_type distance = 0;
    vidType start = rowptr[source];
    dir = Direction::TOP_DOWN;
    if (is_unconnected(source)) {
      distances[source] = 0;
      return;
    } else if (is_leaf(source)) {
      vidType neighbor = copy_unmarked(start);
      if (is_leaf(neighbor)) {
        // Source and only neighbor are both leaf nodes, end BFS
        distances[source] = 0;
        distances[neighbor] = 1;
        return;
      } else {
        this_frontier.push_back(rowptr[neighbor]);
        mark(rowptr[neighbor]);
        distance = 1;
      }
    } else {
      this_frontier.push_back(start);
      mark(start);
    }
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      if (dir == Direction::TOP_DOWN) {
        if (switch_to_bottom_up(&this_frontier)) {
          dir = Direction::BOTTOM_UP;
          bottom_up_step(this_frontier, &next_frontier, distance);
        } else {
          top_down_step(this_frontier, &next_frontier, distance);
        }
      } else {
        if (dir == Direction::BOTTOM_UP) {
          if (switch_to_top_down(&this_frontier)) {
          dir = Direction::TOP_DOWN;
          top_down_step(this_frontier, &next_frontier, distance);
          } else {
            bottom_up_step(this_frontier, &next_frontier, distance);
          }
        }
      }
      // top_down_step(this_frontier, &next_frontier, distance);
      distance++;
      std::swap(this_frontier, next_frontier);
    }
    LIKWID_MARKER_STOP("BFS");
    LIKWID_MARKER_CLOSE;
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "BFS: " << ms_double.count() << "ms\n";

    t1 = std::chrono::high_resolution_clock::now();
    compute_distances(distances, source);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "Postprocessing: " << ms_double.count() << "ms\n";
  }
};

void merged_csr(eidType *rowptr, vidType *col, degrees_map *degrees, uint64_t N, uint64_t M) {
  
  for (vidType i = 0; i < N; i++) {
    // Check if node has more than one neighbor, otherwise don't replace
    // neighbor IDs with indices in col because it makes the computation
    // of distances more complicated
    if (rowptr[i + 1] - rowptr[i] > 1) {
      // Replace list of neighbors with list of neighbor indices in col
      for (uint64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
        col[j] = vidType(rowptr[col[j]]);
      }
    }
    // Mark last neighbor of node
    col[rowptr[i + 1] - 1] = col[rowptr[i + 1] - 1] | MARKED;

    degrees->insert({rowptr[i], rowptr[i + 1] - rowptr[i]});
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  auto t1 = std::chrono::high_resolution_clock::now();
  degrees_map *degrees = new degrees_map();
  degrees->reserve(N);
  merged_csr(rowptr, col, degrees, N, M);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Preprocessing: " << ms_double.count() << "ms\n";
  return new Graph(rowptr, col, degrees, N, M);
}

} // namespace bfs