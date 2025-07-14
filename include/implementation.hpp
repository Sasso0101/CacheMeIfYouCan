#pragma once
#include <cstdint>
#include <vector>
#include "../distributed_mmio/include/mmio.h"

typedef uint32_t vidType;
typedef uint32_t eidType;
typedef uint32_t weight_type;
typedef CSR_local<uint32_t, float> Graph;

typedef enum { TOP_DOWN, BOTTOM_UP } Direction;

using frontier = std::vector<eidType>;

#define ALPHA 4
#define BETA 24

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