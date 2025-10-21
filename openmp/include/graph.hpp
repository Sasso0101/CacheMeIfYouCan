#pragma once
#include <cstdint>
#include <vector>
#include "mmio.h"

typedef uint32_t vertex;
typedef uint32_t edge;

typedef enum { TOP_DOWN, BOTTOM_UP } Direction;

using frontier = std::vector<edge>;

#define ALPHA 4
#define BETA 24

// Graph class to store the graph representation in CSR format
class MergedCSR {
private:
public:
  edge *row_ptr;
  vertex *col_idx;
  uint64_t nrows;
  uint64_t nnz;
};

// Base class for BFS implementations
class BFS_Impl {
public:
  const CSR_local<uint32_t, float> *graph;
  virtual void BFS(vertex source, uint32_t *distances) = 0;
  virtual bool check_result(vertex source, uint32_t *distances) = 0;
  bool check_distances(vertex source, const uint32_t *distances) const;
  bool check_parents(vertex source, const uint32_t *parents) const;
  virtual ~BFS_Impl() = default;

protected:
  BFS_Impl(const CSR_local<uint32_t, float> *graph) : graph(graph) {}
};

// BFS implementation using the MergedCSR graph representation
class MergedCSR_Distances : public BFS_Impl {
private:
  edge *merged_rowptr;
  edge *merged_csr;

  void top_down_step(const frontier &this_frontier, frontier &next_frontier,
                     const uint32_t &distance);
  void compute_distances(uint32_t *distances) const;
  void create_merged_csr();

public:
  MergedCSR_Distances(const CSR_local<uint32_t, float> *graph);
  ~MergedCSR_Distances();
  void BFS(vertex source, uint32_t *distances) override;
  bool check_result(vertex source, uint32_t *distances) override;
};

// BFS implementation using the MergedCSR graph representation (returning
// parents)
class MergedCSR_Parents : public BFS_Impl {
private:
  edge *merged_rowptr;
  edge *merged_csr;

  void top_down_step(const frontier &this_frontier, frontier &next_frontier);
  void compute_parents(uint32_t *parents) const;
  void create_merged_csr();

public:
  MergedCSR_Parents(const CSR_local<uint32_t, float> *graph);
  ~MergedCSR_Parents();
  void BFS(vertex source, uint32_t *distances) override;
  bool check_result(vertex source, uint32_t *distances) override;
};

// Single-threaded BFS implementation using classic CSR
class Reference : public BFS_Impl {
public:
  Reference(const CSR_local<uint32_t, float> *graph);
  ~Reference();
  void BFS(vertex source, uint32_t *distances) override;
  bool check_result(vertex source, uint32_t *distances) override;
};