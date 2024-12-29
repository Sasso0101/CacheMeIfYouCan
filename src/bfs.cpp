#include <chrono>
#include <graph.hpp>
#include <iostream>
#include <omp.h>
#include <parlay/sequence.h>

namespace bfs {

#define MARKED 1 << 31
#define MSB 31

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {}
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

  void BFS(vidType source, weight_type *distances) {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<vidType> this_frontier;
    weight_type distance = 0;
    vidType start = rowptr[source];
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
        distance = 1;
      }
    } else {
      this_frontier.push_back(start);
    }
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto &v : this_frontier) {
        vidType curr_index = v;
        vidType neighbor_start = copy_unmarked(curr_index);
        // Repeat until all neighbors have been visited except last one
        do {
          if (!is_marked(neighbor_start)) {
            mark(neighbor_start);
            next_frontier.push_back(neighbor_start);
          }
          curr_index++;
          neighbor_start = copy_unmarked(curr_index);
        } while (!is_marked(curr_index));
        // Visit last neighbor
        if (!is_marked(neighbor_start)) {
          mark(neighbor_start);
          next_frontier.push_back(neighbor_start);
        }
        set_distance(v, distance);
      }
      distance++;
      std::swap(this_frontier, next_frontier);
    }
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

void merged_csr(eidType *rowptr, vidType *col, uint64_t N, uint64_t M) {
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
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  auto t1 = std::chrono::high_resolution_clock::now();
  merged_csr(rowptr, col, N, M);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Preprocessing: " << ms_double.count() << "ms\n";
  return new Graph(rowptr, col, N, M);
}

} // namespace bfs