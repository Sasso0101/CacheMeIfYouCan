#include "graph.hpp"
#include <cstdint>
#include <gtest/gtest.h>

class BFSTest : public testing::Test {
 protected:
  BFSTest() {
    chdir("../../");
    std::string schema_path = std::string("schemas/Collaboration_Network_1.json");
    g = new Graph(schema_path);
  }
  
  Graph *g;
};

void test_implementation(BFS_Impl *impl, uint32_t source) {
  weight_type *result = new weight_type[impl->graph->N];
  std::fill_n(result, impl->graph->N, std::numeric_limits<weight_type>::max());
  impl->BFS(source, result);
  bool correct = impl->check_result(source, result);
  delete[] result;
  EXPECT_TRUE(correct);
}

TEST_F(BFSTest, Bitmap) {
  BFS_Impl *bitmap = new Bitmap(g);
  test_implementation(bitmap, 5);
}

TEST_F(BFSTest, MergedCSR) {
  BFS_Impl *merged_csr = new MergedCSR(g);
  test_implementation(merged_csr, 5);
}

TEST_F(BFSTest, MergedCSR_1) {
  BFS_Impl *merged_csr_1 = new MergedCSR_1(g);
  test_implementation(merged_csr_1, 5);
}

TEST_F(BFSTest, MergedCSR_Parents) {
  BFS_Impl *mergedCSR_Parents = new MergedCSR_Parents(g);
  test_implementation(mergedCSR_Parents, 5);
}

TEST_F(BFSTest, Classic) {
  BFS_Impl *classic = new Classic(g);
  test_implementation(classic, 5);
}

TEST_F(BFSTest, Reference) {
  BFS_Impl *reference = new Reference(g);
  test_implementation(reference, 5);
}