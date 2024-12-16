#pragma once
#include <cstdint>

typedef uint32_t vidType;
typedef uint64_t eidType;
typedef uint32_t weight_type;

typedef struct {
  vidType parent;
  vidType node;
} pruned;

class BaseGraph {
public:
  virtual ~BaseGraph() {}
  virtual void BFS(vidType source, weight_type *distances) = 0;
};