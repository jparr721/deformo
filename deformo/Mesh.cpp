#include "Mesh.h"

int Vertex::PositionOffset() { return offsetof(Vertex, position); }
int Vertex::ColorOffset() { return offsetof(Vertex, color); }
int Vertex::PositionSize() { return 3; }
int Vertex::ColorSize() { return 3; }
int Vertex::Stride() { return sizeof(Vertex); }

Mesh::Mesh(const Eigen::VectorXd& vertices,
           const std::vector<QVector3D>& colors_)
    : raw_positions(vertices), colors(colors_) {
  LoadVertices();
  IndexDuplicateVertices();
}

void Mesh::LoadVertices() {
  for (auto i = 0u, j = 0u; i < raw_positions.rows();
       i += kNumDimensions, ++j) {
    positions.push_back(Vertex(
        QVector3D(raw_positions(i), raw_positions(i + 1), raw_positions(i + 2)),
        colors[j]));
  }
}

void Mesh::IndexDuplicateVertices() {
  unsigned int node_number = 1;
  std::unordered_map<xyz, unsigned int, xyz_hash> indices;

  for (auto i = 0u; i < raw_positions.rows(); i += kNumDimensions) {
    const auto x = raw_positions(i);      // x
    const auto y = raw_positions(i + 1);  // y
    const auto z = raw_positions(i + 2);  // z
    const xyz pos{x, y, z};

    const bool found = indices.find(pos) != indices.end();

    const unsigned int index = found ? indices.at(pos) : node_number;
    indices.insert({pos, index});

    // Add node number to the vertex position
    positions[i / kNumDimensions].node_number = node_number;

    if (!found) {
      ++node_number;
    }
  }
}

void Mesh::UpdatePositions(Eigen::Ref<const Eigen::VectorXd> displacements) {
  // Nodal displacements are a vector of changes to the nodes after time
  // integration
  raw_positions *= displacements;
  ReloadVertices();
}

void Mesh::ReloadVertices() {
  for (auto i = 0u, j = 0u; i < positions.size(); ++i, j += kNumDimensions) {
    const auto x = raw_positions(j);
    const auto y = raw_positions(j + 1);
    const auto z = raw_positions(j + 2);

    positions[i].position.setX(x);
    positions[i].position.setY(y);
    positions[i].position.setZ(z);
  }
}

void Mesh::Initialize(QOpenGLBuffer& vbo) {
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(positions.data(), positions.size() * sizeof(Vertex));
}

void Mesh::Render(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao) {
  // Write latest vertex positions, we can cache these later if none changed.
  vbo.bind();
  vbo.write(Vertex::PositionOffset(), positions.data(),
            positions.size() * sizeof(Vertex));
  vbo.release();

  vao.bind();
  glDrawArrays(GL_TRIANGLES, 0, positions.size());
  vao.release();
}
