#include <graph.hpp>

namespace test_openmp {
BaseGraph *initialize_graph(eidType *rowptr, vidType *col, uint64_t N,
                            uint64_t M);
}