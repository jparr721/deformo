#include "Loader.h"

#include <fstream>
#include <iostream>

void Loader::LoadTetgen(Eigen::MatrixXd& V, Eigen::MatrixXd& F,
                        const std::string node_file,
                        const std::string ele_file) {
  std::ifstream node_input;
  std::ifstream ele_input;

  node_input.open(node_file);
  ele_input.open(ele_file);

  if (!node_input.good()) {
    std::cerr << "Cannot read from " << node_file << std::endl;
    exit(1);
  }

  if (!ele_input.good()) {
    std::cerr << "Cannot read from " << node_file << std::endl;
    exit(1);
  }

  ReadTetgenVertexFile(V, node_input);
  ReadTetgenEleFile(F, ele_input);
}

void Loader::ReadTetgenVertexFile(Eigen::MatrixXd& V,
                                  std::ifstream& node_input) {
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

void Loader::ReadTetgenEleFile(Eigen::MatrixXd& F, std::ifstream& ele_input) {
  int n_nodes = 0;
  int n_vertices = 0;
  std::string line;

  ele_input >> n_nodes;
  F.resize(n_nodes, 4);

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

    // Discard node number
    ss >> node_number;

    // Then iterate the rest of the line
    for (int i = 0; i < n_vertices; ++i) {
      // F (row, column) = value from ss (automatically casted)
      ss >> F(row, i);
    }

    ++row;

    // End of the matrix
    if (row >= F.rows()) {
      break; 
    }
  }
}
