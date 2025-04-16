#pragma once
#include <atomic>
#include <cstdint>
#include <omp.h>
#include <semaphore.hpp>
#include <string>
#include <vector>

typedef uint32_t vidType;
typedef uint32_t eidType;
typedef uint32_t weight_type;

typedef enum {
  merged_csr_parents,
  merged_csr,
  merged_csr_1,
  bitmap,
  classic,
  reference,
  heuristic
} Implementation;

typedef enum { TOP_DOWN, BOTTOM_UP } Direction;

using frontier = std::vector<eidType>;

#define ALPHA 4
#define BETA 24

// Graph class to store the graph representation in CSR format
class Graph {
private:
  void construct_from_coo(std::vector<int64_t> &input_row,
                          std::vector<int64_t> &input_col);
  void construct_from_file(std::string &filename);
  void generate_random_graph(int64_t num_vertices,
                             int64_t num_edges_per_vertex);

public:
  eidType *rowptr;
  vidType *col;
  uint64_t N;
  uint64_t M;

  Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  Graph(std::string &filename);
  ~Graph();
  void print_graph();
};

// Base class for BFS implementations
class BFS_Impl {
public:
  Graph *graph;
  virtual void BFS(vidType source, weight_type *distances) = 0;
  virtual bool check_result(vidType source, weight_type *distances) = 0;
  bool check_distances(vidType source, const weight_type *distances) const;
  bool check_parents(vidType source, const weight_type *parents) const;

protected:
  BFS_Impl(Graph *graph) : graph(graph) {}
  ~BFS_Impl() { delete graph; }
};

// BFS implementation using bitmaps to store frontiers and visited array
class Bitmap : public BFS_Impl {
private:
  bool *this_frontier;
  bool *next_frontier;
  bool *visited;

  void bottom_up_step(const bool *this_frontier, bool *next_frontier);
  void top_down_step(const bool *this_frontier, bool *next_frontier);
  inline void add_to_frontier(bool *frontier, vidType v);

public:
  Bitmap(Graph *graph);
  ~Bitmap();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};

// BFS implementation using the MergedCSR graph representation
class MergedCSR : public BFS_Impl {
private:
  eidType *merged_rowptr;
  eidType *merged_csr;

  void top_down_step(const frontier &this_frontier, frontier &next_frontier,
                     const weight_type &distance);
  void compute_distances(weight_type *distances, vidType source) const;
  void create_merged_csr();

public:
  MergedCSR(Graph *graph);
  ~MergedCSR();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};

// BFS implementation using the MergedCSR graph representation (returning
// parents)
class MergedCSR_Parents : public BFS_Impl {
private:
  eidType *merged_rowptr;
  eidType *merged_csr;

  void top_down_step(const frontier &this_frontier, frontier &next_frontier);
  void compute_parents(weight_type *parents, vidType source) const;
  void create_merged_csr();

public:
  MergedCSR_Parents(Graph *graph);
  ~MergedCSR_Parents();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};

// BFS implementation using bitmaps to store visited array. Frontiers are stored
// as vectors
class Classic : public BFS_Impl {
private:
  bool *visited;

  inline void set_distance(vidType i, weight_type distance,
                           weight_type *distances);
  inline void add_to_frontier(frontier &frontier, vidType v,
                              vidType &edges_frontier);
  void bottom_up_step(frontier this_frontier, frontier &next_frontier,
                      weight_type distance, weight_type *distances,
                      vidType &edges_frontier);
  void top_down_step(frontier this_frontier, frontier &next_frontier,
                     weight_type &distance, weight_type *distances,
                     vidType &edges_frontier, vidType edges_frontier_old);

public:
  Classic(Graph *graph);
  ~Classic();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};

// Single-threaded BFS implementation using classic CSR
class Reference : public BFS_Impl {
public:
  Reference(Graph *graph);
  ~Reference();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};

// BFS implementation using the MergedCSR graph representation
class MergedCSR_1 : public BFS_Impl {
private:
  eidType *merged_rowptr;
  eidType *merged_csr;
  std::atomic<std::uint32_t> waiting_threads;
  std::vector<vidType> comm_frontier;
  weight_type comm_distance;
  semaphore *readers_lock;
  semaphore *comm_lock;
  uint64_t overwrites;
  uint64_t reduntants;
  std::vector<uint32_t> thread_errors;

  void top_down_step(const frontier &this_frontier, frontier &next_frontier,
                     const weight_type &distance);
  void compute_distances(weight_type *distances, vidType source) const;
  void parallel_frontiers(const frontier &this_frontier, weight_type &distance);
  void create_merged_csr();
  void check_vertex(eidType i, frontier &next_frontier,
                    const weight_type &local_distance);
  void critical_writeback(frontier &next_frontier,
                          const weight_type &local_distance);
  void balance_threads(frontier &next_frontier, weight_type &local_distance);
  void donate_frontier(frontier &next_frontier, weight_type &local_distance);
  void obtain_frontier(frontier &next_frontier, weight_type &local_distance);

public:
  MergedCSR_1(Graph *graph);
  ~MergedCSR_1();
  void BFS(vidType source, weight_type *distances) override;
  bool check_result(vidType source, weight_type *distances) override;
};