#include <graph.hpp>

namespace benchmark {
BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M);
}