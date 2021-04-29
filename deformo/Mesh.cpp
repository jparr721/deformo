#include "Mesh.h"

#include <iostream>
Mesh::Mesh(const Eigen::VectorXd& vertices_) : vertices(vertices_) {
  LoadVBO();
}

void Mesh::LoadVBO() {
  unsigned int node_number = 1;
  for (auto i = 0u; i < vertices.rows(); i += 2) {
    const xy pos{
        vertices(i),      // x
        vertices(i + 1),  // y
    };

    const bool found = indices.find(pos) != indices.end();

    const unsigned int index = found ? indices.at(pos) : node_number;
    indices.insert({{pos, index}});

    if (!found) {
      ++node_number;
    }
  }
}
