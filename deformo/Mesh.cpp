#include "Mesh.h"

Mesh::Mesh(const Eigen::VectorXd& vertices_) : vertices(vertices_) {
  LoadVBO();
}

void Mesh::LoadVBO() {
  unsigned int node_number = 0;
  for (auto i = 0u; i < vertices.rows(); i += 2) {
    const xy pos{
        vertices(i),      // x
        vertices(i + 1),  // y
    };

    const unsigned int index =
        indices.find(pos) != indices.end() ? indices.at(pos) : node_number;

    indices.insert({{pos, index}});
    ++node_number;
  }
}
