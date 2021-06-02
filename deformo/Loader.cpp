#include "Loader.h"

#include <fstream>

namespace loader {

namespace {
std::ifstream OpenFile(const std::string& filename) {
  std::ifstream input;
  input.open(filename);
  assert(input.good());

  return input;
}
}  // namespace

void ReadTetgenVertexFile(Eigen::MatrixXf& V, const std::string& node_file) {
  auto node_input = OpenFile(node_file);
  int n_nodes;
  node_input >> n_nodes;

  // Resize to fit the number of nodes in the vertex file
  V.resize(n_nodes, 3);

  unsigned int row = 0;

  std::string line;

  // Get rid of the rest of the first line
  std::getline(node_input, line);

  int node_number = 0;

  for (;;) {
    std::getline(node_input, line);

    if (node_input.eof()) {
      break;
    }

    if (line.size() < 3) {
      continue;
    }

    // Skip Comments
    if (line.at(0) == '#') {
      continue;
    }

    std::stringstream ss(line);

    // Collect the node number first
    ss >> node_number;

    // Then iterate the rest of the line
    for (int i = 0; i < 3; ++i) {
      // V (row, column) = value from ss (automatically casted)
      ss >> V(row, i);
    }

    ++row;

    // End of the matrix
    if (row >= V.rows()) {
      break;
    }
  }
}

void ReadTetgenFaceFile(Eigen::MatrixXi& F, const std::string& face_file) {
  auto face_input = OpenFile(face_file);
  int n_nodes = 0;
  int row = 0;
  std::string line;

  face_input >> n_nodes;

  // Ignore boundary faces for now
  F.resize(n_nodes, 3);

  // Discard the rest of the line
  std::getline(face_input, line);

  for (;;) {
    std::getline(face_input, line);

    // Tetgen has the last line as a comment
    if (face_input.eof()) {
      break;
    }

    if (line.size() < 3) {
      continue;
    }

    // Skip comments
    if (line.at(0) == '#') {
      continue;
    }

    std::stringstream ss(line);

    int _buf_node_number;
    // Discard node number, index should keep track of it just fine
    ss >> _buf_node_number;

    for (int i = 0; i < 3; ++i) {
      ss >> F(row, i);
    }

    ++row;

    if (row >= F.rows()) {
      break;
    }
  }
}

void ReadTetgenEleFile(Eigen::MatrixXi& T, const std::string& ele_file) {
  auto ele_input = OpenFile(ele_file);
  int n_nodes = 0;
  int n_vertices = 0;
  std::string line;

  ele_input >> n_nodes;
  T.resize(n_nodes, 4);

  ele_input >> n_vertices;

  unsigned int row = 0;
  int node_number = 0;

  // Discard the rest of the line
  std::getline(ele_input, line);

  for (;;) {
    std::getline(ele_input, line);

    if (ele_input.eof()) {
      break;
    }

    if (line.size() < 3) {
      continue;
    }

    // Skip Comments
    if (line.at(0) == '#') {
      continue;
    }

    std::stringstream ss(line);

    // Discard node number, index should keep track of it just fine
    ss >> node_number;

    // Then iterate the rest of the line
    for (int i = 0; i < n_vertices; ++i) {
      // T (row, column) = value from ss (automatically casted)
      ss >> T(row, i);
    }

    ++row;

    // End of the matrix
    if (row >= T.rows()) {
      break;
    }
  }
}
}  // namespace loader
