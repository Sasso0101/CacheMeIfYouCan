#include "graph.hpp"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <random>
#include <string>
#include "inputschema.cpp"

Graph::Graph(eidType *rowptr, vidType *col, uint64_t N, uint64_t M)
    : rowptr(rowptr), col(col), N(N), M(M) {}

Graph::Graph(std::string &filename) {
  nlohmann::json j;
  std::ifstream in(filename);
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  if (data.graph.data_file_format.has_value()) {
    assert(data.graph.filename.has_value() &&
            data.graph.file_format.has_value());
    assert(data.graph.file_format.value() ==
            "binary");
    construct_from_file(data.graph.filename.value());
  } else if (data.graph.coo_format.has_value()) {
    assert(data.graph.row.has_value() && data.graph.col.has_value()
            && "COO values missing.");
    construct_from_coo(data.graph.row.value(), data.graph.col.value());
  } else if (data.graph.random_generated_graph.has_value()) {
    generate_random_graph(data.graph.num_vertices.value(),
                          data.graph.num_edges_per_vertex.value());
  } else {
    assert(false && "Error no valid format\n");
  }
}

Graph::~Graph() {
  delete[] rowptr;
  delete[] col;
}

void Graph::construct_from_coo(std::vector<int64_t> &input_row,
                                   std::vector<int64_t> &input_col) {
  // convert COO form to CSR format.
  N = input_row[0];
  for (int64_t i = input_row.size(); --i >= 0;) {
    if (input_row[i] + 1 > N) {
      N = input_row[i] + 1;
    }
  }
  M = input_col.size();
  assert(input_col.size() == input_row.size() &&
         "In COO format col and row must have the same lengths");
  rowptr = new eidType[N + 1]();
  col = new vidType[M]();

  for (size_t i = 0; i < input_row.size(); i++) {
    vidType src = input_row[i];
    vidType dst = input_col[i];
    col[i] = dst;
    rowptr[src + 1] = i + 1;
  }
}

void Graph::construct_from_file(std::string &filename) {
  std::string path = "datasets/" + filename;
  std::ifstream s{path, s.in | s.binary};
  if (!s.is_open()) {
    throw std::runtime_error("Error: Unable to open file " + path);
  }

  s.read((char *)&N, sizeof(decltype(N)));
  s.read((char *)&M, sizeof(decltype(M)));

  rowptr = new eidType[N + 1];
  col = new vidType[M];

  uint64_t *temp_rowptr = new uint64_t[N + 1];
  s.read((char *)temp_rowptr, sizeof(uint64_t) * (N + 1));
  // Convert to eidType from uint64_t
  for (uint64_t i = 0; i <= N; i++) {
    rowptr[i] = static_cast<eidType>(temp_rowptr[i]);
  }
  delete[] temp_rowptr;
  s.read((char *)col, sizeof(uint32_t) * M);
  
  s.close();
}

void Graph::generate_random_graph(int64_t num_vertices,
                                      int64_t num_edges_per_vertex) {
  std::random_device r;
  std::default_random_engine el(r());
  std::uniform_int_distribution<uint64_t> uniform_dist(0, num_vertices);
  int64_t num_edges = num_vertices * num_edges_per_vertex;
  std::vector<std::tuple<int64_t, int64_t>> edges;
  edges.reserve(num_edges);
  for (int64_t i = 0; i < num_edges; i++) {
    uint64_t src = uniform_dist(el);
    uint64_t dst = uniform_dist(el);
    if (src == dst)
      continue;
    edges.push_back(std::make_tuple(src, dst));
    edges.push_back(std::make_tuple(dst, src));
  }
  std::sort(edges.begin(), edges.end());
  std::vector<std::tuple<int64_t, int64_t>> filtered_edges;
  filtered_edges.reserve(edges.size());
  for (int64_t i = 0; i < edges.size(); i++) {
    if (i == 0 || std::get<0>(edges[i]) != std::get<0>(edges[i - 1]) ||
        std::get<1>(edges[i]) != std::get<1>(edges[i - 1])) {
      filtered_edges.push_back(edges[i]);
    }
  }
  std::vector<int64_t> row(filtered_edges.size());
  std::vector<int64_t> col(filtered_edges.size());
  for (int64_t i = 0; i < filtered_edges.size(); i++) {
    row[i] = std::get<0>(filtered_edges[i]);
    col[i] = std::get<1>(filtered_edges[i]);
  }
  construct_from_coo(row, col);
}

void Graph::print_graph() {
  // Print the number of vertices (N) and edges (M) in the graph
  printf("Number of vertices (N): %lu\n", N);
  printf("Number of edges (M): %lu\n", M);
  // Print the rowptr array of the graph
  printf("Rowptr: ");
  for (size_t i = 0; i < 10; ++i) {
    printf("%u ", rowptr[i]);
  }
  printf("\n");

  // Print the colidx array of the graph
  printf("Col: ");
  for (size_t i = 0; i < 10; ++i) {
    printf("%u ", col[i]);
  }
  printf("\n");
}