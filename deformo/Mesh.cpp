#include "Mesh.h"

int Vertex::PositionOffset() { return offsetof(Vertex, position); }
int Vertex::ColorOffset() { return offsetof(Vertex, color); }
int Vertex::PositionSize() { return 3; }
int Vertex::ColorSize() { return 3; }
int Vertex::Stride() { return sizeof(Vertex); }

Mesh::Mesh(const Eigen::VectorXd& vertices, bool is_2d)
    : raw_positions(vertices) {
  LoadVertices(vertices, is_2d);
  IndexDuplicateVertices();
}

Mesh::Mesh(const Eigen::VectorXd& vertices,
           const std::vector<QVector3D>& colors)
    : raw_positions(vertices) {
  LoadVertices(vertices, colors);
  IndexDuplicateVertices();
}

void Mesh::LoadVertices(const Eigen::VectorXd& vertices, bool is_2d) {
  if (is_2d) {
    for (int i = 0u; i < raw_positions.rows(); i += 2) {
      positions.push_back(
          Vertex(QVector3D(raw_positions(i), raw_positions(i + 1), -3.)));
    }
  } else {
    for (auto i = 0u; i < raw_positions.rows(); i += 3) {
      positions.push_back(Vertex(QVector3D(
          raw_positions(i), raw_positions(i + 1), raw_positions(i + 2))));
    }
  }
}

void Mesh::LoadVertices(const Eigen::VectorXd& vertices,
                        const std::vector<QVector3D>& colors) {
  for (auto i = 0u, j = 0u; i < raw_positions.rows(); i += 3, ++j) {
    positions.push_back(Vertex(
        QVector3D(raw_positions(i), raw_positions(i + 1), raw_positions(i + 2)),
        colors[j]));
  }
}

void Mesh::IndexDuplicateVertices() {
  unsigned int node_number = 1;
  for (auto i = 0u; i < raw_positions.rows(); i += 2) {
    const auto x = raw_positions(i);      // x
    const auto y = raw_positions(i + 1);  // y
    const xy pos{x, y};

    const bool found = indices.find(pos) != indices.end();

    const unsigned int index = found ? indices.at(pos) : node_number;
    indices.insert({{pos, index}});

    if (!found) {
      ++node_number;
    }
  }
}

void Mesh::Initialize(QOpenGLBuffer& vbo) {
  vbo.create();
  vbo.bind();
  vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
  vbo.allocate(positions.data(), positions.size() * sizeof(Vertex));
}

void Mesh::Render(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao) {
  vbo.bind();
  vbo.write(Vertex::PositionOffset(), positions.data(),
            positions.size() * sizeof(Vertex));
  vbo.release();
  vao.bind();
  glDrawArrays(GL_TRIANGLES, 0, positions.size());
  vao.release();
}
