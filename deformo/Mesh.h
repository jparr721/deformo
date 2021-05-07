#pragma once

#include <Eigen/Dense>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QVector3D>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

constexpr unsigned int kTetrahedronElementCount = 4;

struct xyz {
  double x;
  double y;
  double z;
  bool operator==(const xyz& rhs) const {
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
};

struct xyz_hash {
  std::size_t operator()(const xyz& self) const noexcept {
    std::size_t seed = 3;

    seed ^= static_cast<unsigned int>(self.x) + 0x9e779b9 + (seed << 6) +
            (seed >> 2);
    seed ^= static_cast<unsigned int>(self.y) + 0x9e779b9 + (seed << 6) +
            (seed >> 2);
    seed ^= static_cast<unsigned int>(self.z) + 0x9e779b9 + (seed << 6) +
            (seed >> 2);

    return seed;
  }
};

struct Vertex {
  QVector3D position;
  QVector3D color;

  Vertex(const QVector3D& vertex_)
      : position(vertex_), color(QVector3D(1., 0., 0.)) {}
  Vertex(const QVector3D& vertex_, const QVector3D& color_)
      : position(vertex_), color(color_) {}

  static int PositionOffset();
  static int ColorOffset();
  static int PositionSize();
  static int ColorSize();
  static int Stride();
};

class Mesh {
 public:
  static constexpr unsigned int kNumDimensions = 3;

  std::vector<Vertex> positions;
  std::vector<QVector3D> colors;

  std::unordered_map<xyz, unsigned int, xyz_hash> indices;

  Eigen::VectorXd raw_positions;

  Mesh(const Eigen::VectorXd& vertices, const std::vector<QVector3D>& colors_);

  void UpdatePositions(Eigen::Ref<const Eigen::VectorXd> displacements);
  void Initialize(QOpenGLBuffer& vbo);
  void Render(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao);

  unsigned int rows() { return raw_positions.rows(); }
  unsigned int index(double x, double y, double z) {
    return indices.at(xyz{x, y, z});
  }

 private:
  void LoadVertices();
  void ReloadVertices();

  void IndexDuplicateVertices();
};
