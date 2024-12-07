#include <bitset>
#include <cmath>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>

typedef uint64_t vidType;
typedef uint64_t eidType;
typedef uint64_t weight_type;

class BaseGraph {
public:
  virtual ~BaseGraph() {}
  virtual void BFS(vidType source, weight_type *distances) = 0;
};

class Graph : public BaseGraph {
  eidType *rowptr;
  vidType *col;
  [[maybe_unused]] uint64_t N;
  [[maybe_unused]] uint64_t M;

public:
  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M) {}
  ~Graph() {
    // destructor logic.
    // If you perform any memory allocations with malloc, new, etc. you must
    // free
    //   them here to avoid memory leaks.
  }

  void BFS(vidType source, weight_type *distances) {
    std::vector<vidType> this_frontier;
    distances[source] = 0;
    this_frontier.push_back(source);
    while (!this_frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto &src : this_frontier) {
        for (uint64_t i = rowptr[src]; i < rowptr[src + 1]; i++) {
          vidType dst = col[i];
          if (distances[src] + 1 < distances[dst]) {
            distances[dst] = distances[src] + 1;
            next_frontier.push_back(dst);
          }
        }
      }
      std::swap(this_frontier, next_frontier);
    }
  }
};

template <typename T>
uint64_t binary_search(T *arr, uint64_t start, uint64_t end, T value) {
  uint64_t mid = (start + end) / 2;
  while (start < end) {
    if (arr[mid] == value) {
      return mid;
    } else if (arr[mid] < value) {
      start = mid + 1;
    } else {
      end = mid;
    }
    mid = (start + end) / 2;
  }
  return start;
}

void remove_node_degree_1(eidType *rowptr, vidType *col, vidType v, eidType neighbor) {
  // Remove the node from the neighbor's adjacency list
  uint64_t to_remove = binary_search(col, rowptr[neighbor], rowptr[neighbor + 1]-1, v);
  std::cout << "Start " << rowptr[neighbor] << " end " << rowptr[neighbor + 1]-1 << std::endl;
  std::cout << "Removed node " << v << " marking " << col[to_remove] << std::endl;
  col[to_remove] = -1;
}

void print_graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M) {
  for (uint64_t i = 0; i < N; i++) {
    for (uint64_t j = rowptr[i]; j < rowptr[i+1]; j++) {
      std::cout << i << " (j: " << j << ") " << col[j] << std::endl;
    }
  }
}

BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M) {
  std::unordered_multiset<eidType> to_check;
  uint64_t removed = 0;
  // Remove nodes with only one neighbor
  for (uint64_t v = 0; v < N; v++) {
    if (rowptr[v + 1] - rowptr[v] == 1) {
      eidType neighbor = col[rowptr[v]];
      // Insert the neighbor into the set of nodes because it may be removed
      to_check.insert(neighbor);
      remove_node_degree_1(rowptr, col, v, neighbor);
      removed++;
    }
  }

  // Continue removing nodes with only one neighbor
  while (to_check.size() > 0) {
    vidType v = *to_check.begin();
    uint64_t removed_neighbors = to_check.count(v);
    to_check.erase(v);
    if (rowptr[v+1] - rowptr[v] == removed_neighbors + 1) {
      removed++;
      // Find the neighbor
      for (uint64_t i = rowptr[v]; i < rowptr[v + 1]; i++) {
        if (col[i] != -1) {
          vidType neighbor = col[i];
          to_check.insert(neighbor);
          remove_node_degree_1(rowptr, col, v, neighbor);
          break;
        }
      }
    }
  }
  std::cout << "Removed " << removed << " nodes" << std::endl;

  return new Graph(rowptr, col, N, M);
}