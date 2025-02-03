#include "../inputschema.cpp"
#include <fstream>
#include <graph.hpp>
#include <string>

#define SEED 458509012

class ProblemInput {
  quicktype::Inputschema input;
  uint64_t N;
  uint64_t M;

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
  }

  void construct_from_file(std::string &filename) {
    // std::ifstream s{filename, s.trunc | s.in | s.out | s.binary};
    std::string path = "datasets/" + filename;
    std::ifstream s{path, s.in | s.binary};

    s.read((char *)&N, sizeof(decltype(N)));
    s.read((char *)&M, sizeof(decltype(M)));
  }

  ProblemInput(quicktype::Inputschema &_input, int num_sources, std::string name) {
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
    } else {
      assert(false && "Error no valid format\n");
    }

    // Open sources.txt by appending to it at the end
    std::ofstream source_file("sources.txt", std::ios_base::app);
    srand(SEED);
    source_file << name << " ";
    for (int i = 0; i < num_sources; i++) {
      source_file << rand() % N << " ";
    }
    source_file << std::endl;
  }
};

#define USAGE "Usage: %s <schema> num_sources\nThe second argument is the number of sources to be generated.\n"

int main(const int argc, char **argv) {
  if (argc != 3) {
    printf(USAGE, argv[0]);
    return 1;
  }
  std::ifstream in("schemas/" + std::string(argv[1]) + ".json");
  nlohmann::json j;
  in >> j;
  quicktype::Inputschema data;
  quicktype::from_json(j, data);
  ProblemInput p = ProblemInput(data, std::stoi(argv[2]), argv[1]);
}