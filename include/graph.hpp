#pragma once
#include <cstdint>
#include <vector>

typedef uint32_t vidType;
typedef uint32_t eidType;
typedef uint32_t weight_type;

typedef enum {
  BITMAP,
  MERGED_CSR,
  MERGED_CSR_PARENTS,
  CLASSIC,
  REFERENCE,
  HEURISTIC
} Algorithm;
typedef enum { TOP_DOWN, BOTTOM_UP } Direction;

using frontier = std::vector<eidType>;

#define ALPHA 4
#define BETA 24

// Base class for graph representations using classic CSR
class BaseGraph {
protected:
  eidType *rowptr;
  vidType *col;
  uint64_t N;
  uint64_t M;

public:
  BaseGraph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
      : rowptr(rowptr), col(col), N(N), M(M){};
  virtual void BFS(vidType source, weight_type *distances) = 0;
};

// BFS implementation using bitmaps to store frontiers and visited array
class Bitmap : public BaseGraph {
private:
  bool *this_frontier;
  bool *next_frontier;
  bool *visited;

  void bottom_up_step(const bool *this_frontier, bool *next_frontier);
  void top_down_step(const bool *this_frontier, bool *next_frontier);
  inline void add_to_frontier(bool *frontier, vidType v);

public:
  Bitmap(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  void BFS(vidType source, weight_type *distances) override;
};

// BFS implementation using the MergedCSR graph representation
class MergedCSR : public BaseGraph {
private:
  eidType *merged;

  inline void add_to_frontier(frontier &frontier, eidType v) const;
  void top_down_step(const frontier &this_frontier, frontier &next_frontier,
                     const weight_type &distance);
  void compute_distances(weight_type *distances, vidType source) const;
  void create_merged_csr(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);

public:
  MergedCSR(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  void BFS(vidType source, weight_type *distances) override;
};

// BFS implementation using the MergedCSR graph representation (returning
// parents)
class MergedCSR_Parents : public BaseGraph {
private:
  eidType *merged;

  inline void add_to_frontier(frontier &frontier, eidType v) const;
  void top_down_step(const frontier &this_frontier, frontier &next_frontier);
  void compute_parents(weight_type *parents, vidType source) const;
  void create_merged_csr(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);

public:
  MergedCSR_Parents(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  void BFS(vidType source, weight_type *distances) override;
};

// BFS implementation using bitmaps to store visited array. Frontiers are stored
// as vectors
class Classic : public BaseGraph {
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
  Classic(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  void BFS(vidType source, weight_type *distances) override;
};

// Single-threaded BFS implementation using classic CSR
class Reference : public BaseGraph {
public:
  Reference(eidType *rowptr, vidType *col, uint64_t N, uint64_t M);
  void BFS(vidType source, weight_type *distances) override;
};