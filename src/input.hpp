#include "inputschema.cpp"
#include <cassert>
#include <cstdint>
#include <fstream>
#include <graph.hpp>
#include <iostream>
#include <optional>
#include <random>
#include <reference.hpp>
#include <string>
#include <complete.hpp>
#include <parents.hpp>

enum class Problem { Distances, Parents };

class ProblemInput {
  quicktype::Inputschema input;
  uint64_t N;
  uint64_t M;
  eidType *rowptr;
  vidType *col;
  Problem problem;
  std::vector<weight_type> output;
  BaseGraph *graph;

public:
  uint64_t get_num_operations() { return M; }
  void construct_from_coo(std::vector<int64_t> &input_row,
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
      rowptr[src + 1] = i+1;
    }

    std::cout << "N: " << N << " M: " << M << std::endl;
    for (int64_t i = 0; i < N; i++) {
      std::cout << rowptr[i] << " ";
      for (int64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
        std::cout << col[j] << " ";
      }
      std::cout << std::endl;
    }
  }

  void construct_from_file(std::string &filename) {
    // std::ifstream s{filename, s.trunc | s.in | s.out | s.binary};
    std::string path = "datasets/" + filename;
    std::ifstream s{path, s.in | s.binary};

    s.read((char *)&N, sizeof(decltype(N)));
    s.read((char *)&M, sizeof(decltype(M)));

    rowptr =
        new eidType[N + 1]; //(decltype(rowptr)) malloc(sizeof(eidType)*(N+1));
    col = new vidType[M];   //(decltype(col)) malloc(sizeof(vidType)*M);
    // weights = new weight_type[M];//(decltype(weights))
    // malloc(sizeof(weight_type)*M);

    s.read((char *)rowptr, sizeof(eidType) * (N + 1));
    s.read((char *)col, sizeof(vidType) * M);
    // s.read((char*)weights, sizeof(weight_type)*M);
    s.close();
    // printf("Successfully deserialized the data. %llu, %llu\n", N, M);
  }

  ProblemInput(quicktype::Inputschema &_input, std::string algorithm, std::string _problem = "distances") {
    this->input = _input;
    if (_input.graph.data_file_format.has_value()) { // filename.has_value()) {
      assert(_input.graph.filename.has_value() &&
             _input.graph.file_format.has_value());
      assert(_input.graph.file_format.value() ==
             "binary"); // update if there are more file formats.
      construct_from_file(_input.graph.filename.value());
    } else if (_input.graph.coo_format.has_value()) {
      assert(_input.graph.row.has_value() && _input.graph.col.has_value() /*&&
             _input.graph.weight.has_value()*/
             && "COO values missing.");
      construct_from_coo(_input.graph.row.value(), _input.graph.col.value());
    } else if (_input.graph.random_generated_graph.has_value()) {
      generate_random_graph(_input.graph.num_vertices.value(),
                            _input.graph.num_edges_per_vertex.value());
    } else {
      assert(false && "Error no valid format\n");
    }

    // Helper function to convert formats to binary format.
    // NOTE: The serialize function will not replace existing files with the
    // same name in the data/ directory.
    if (_input.meta_info.has_value()) {
      auto meta_info = _input.meta_info.value();
      if (meta_info.save_to_binary) {
        serialize(meta_info.save_filename);
      }
    }

    if (_problem == "distances") {
      this->problem = Problem::Distances;
      graph = complete::initialize_graph(rowptr, col, N, M, algorithm);
      output = std::vector<vidType>(N, std::numeric_limits<weight_type>::max());
    } else if (_problem == "parents") {
      this->problem = Problem::Parents;
      output = std::vector<vidType>(N, -1);
      graph = parents::initialize_graph(rowptr, col, N, M, algorithm);
    } else {
      assert(false && "Invalid problem type");
    }
  }

  void generate_random_graph(int64_t num_vertices,
                             int64_t num_edges_per_vertex) {
    std::random_device r;
    std::default_random_engine el(r());
    std::uniform_int_distribution<uint64_t> uniform_dist(0, num_vertices);
    // std::uniform_real_distribution<float> uniform_weight_dist(0,10.0);
    int64_t num_edges = num_vertices * num_edges_per_vertex;
    std::vector<std::tuple<int64_t, int64_t>> edges;
    edges.reserve(num_edges);
    for (int64_t i = 0; i < num_edges; i++) {
      uint64_t src = uniform_dist(el);
      uint64_t dst = uniform_dist(el);
      if (src == dst)
        continue;
      // float weight = uniform_weight_dist(el);
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
    // std::vector<double> weights(filtered_edges.size());
    for (int64_t i = 0; i < filtered_edges.size(); i++) {
      row[i] = std::get<0>(filtered_edges[i]);
      col[i] = std::get<1>(filtered_edges[i]);
      // weights[i] = std::get<2>(filtered_edges[i]);
    }
    construct_from_coo(row, col);
  }

  void serialize(std::string name) {
    std::string filename = "datasets/" + name; //"test.bin";
    if (std::ifstream(filename.c_str()).good()) {
      printf("Data file already exists, skipping\n");
      return;
    }
    std::ofstream s{filename, s.trunc | s.in | s.out | s.binary};
    s.write((char *)&N, sizeof(decltype(N)));
    s.write((char *)&M, sizeof(decltype(M)));
    s.write((char *)rowptr, sizeof(eidType) * (N + 1));
    s.write((char *)col, sizeof(vidType) * M);
    // s.write((char*)weights, sizeof(weight_type)*M);
    s.close();
  }

  ~ProblemInput() {
    if (rowptr != NULL) {
      delete[] rowptr;
      rowptr = NULL;
    }
    if (col != NULL) {
      delete[] col;
      col = NULL;
    }
    // if (weights != NULL) {
    //	delete[] weights;
    //	weights = NULL;
    // }
    if (graph != NULL) {
      delete graph;
      graph = NULL;
    }
  }

  void release_memory_postrun() {
    if (rowptr != NULL) {
      delete[] rowptr;
      rowptr = NULL;
    }
    if (col != NULL) {
      delete[] col;
      col = NULL;
    }
    // if (weights != NULL) {
    //	delete[] weights;
    //	weights = NULL;
    // }
    if (graph != NULL) {
      delete graph;
      graph = NULL;
    }
  }

  void run(bool check = false, std::optional<vidType> _source = std::nullopt) {
    vidType source;
    if (_source.has_value()) {
      source = _source.value();
    } else {
      source = rand() % N;
    }
    graph->BFS(source, output.data());
    if (check) {
      if (problem == Problem::Distances) {
        std::vector<weight_type> ref_distances(N, std::numeric_limits<weight_type>::max());
        auto ref_graph = reference::initialize_graph(rowptr, col, N, M, "");
        ref_graph->BFS(source, ref_distances.data());
        for (int64_t i = 0; i < output.size(); i++) {
          if (output.size() != ref_distances.size())
            std::cout << "Incorrect # queries in distance array\n";
          for (int64_t j = 0; j < output.size(); j++) {
            if (output[j] != ref_distances[j]) {
              std::cout << "Incorrect value, expected distance " +
                               std::to_string(ref_distances[j]) +
                               ", but got " + std::to_string(output[j]) + "\n";
            }
          }
        }
      } else {
        std::vector<vidType> depth(N, -1);
        std::vector<vidType> to_visit;
        depth[source] = 0;
        to_visit.push_back(source);
        to_visit.reserve(N);
        // Run BFS to compute depth of each vertex
        for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
          vidType i = *it;
          std::cout << "Visiting " << i << std::endl;
          std::cout << "Start " << rowptr[i] << " End " << rowptr[i + 1] << std::endl;
          for (int64_t v = rowptr[i]; v < rowptr[i + 1]; v++) {
            if (depth[col[v]] == -1) {
              depth[col[v]] = depth[i] + 1;
              to_visit.push_back(col[v]);
              std::cout << "Pushing " << col[v] << " with depth " << depth[col[v]] << std::endl;
            }
          }
        }
        // print depth array
        for (int64_t i = 0; i < N; i++) {
          std::cout << output[i] << " ";
        }
        for (int64_t i = 0; i < N; i++) {
          // Check if vertex is part of the BFS tree
          if (depth[i] != -1 && output[i] != -1) {
            // Check if parent is correct
            if (i == source) {
              if (!((output[i] == i) && (depth[i] == 0))) {
                std::cout << "Source wrong";
              }
              continue;
            }
            bool parent_found = false;
            for (int64_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
              if (col[j] == output[i]) {
                vidType parent = col[j];
                // Check if parent has correct depth
                if (depth[parent] != depth[i] - 1) {
                  std::cout << "Wrong depth of child " + std::to_string(i) +
                                   " (parent " + std::to_string(parent) +
                                   "with depth " + std::to_string(depth[parent])
                            << ")" << std::endl;
                }
                parent_found = true;
                break;
              }
            }
            if (!parent_found) {
              std::cout << "Couldn't find edge from " << output[i] << " to "
                        << std::to_string(i) << std::endl;
            }
            // Check if parent = -1 and parent = -1
          } else if (depth[i] != output[i]) {
            std::cout << "Reachability mismatch" << std::endl;
          }
        }
      }
    }
  }
};
