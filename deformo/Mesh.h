#pragma once

#include <Eigen/Dense>
#include <QOpenGLBuffer>
#include <QVector3D>
#include <QOpenGLVertexArrayObject>
#include <unordered_map>
#include <utility>
#include <vector>

struct xy {
  double x;
  double y;
  bool operator==(const xy& other) const {
    return x == other.x && y == other.y;
  }
};

struct xy_hash {
  std::size_t operator()(const xy& instance) const noexcept {
    std::size_t x_hash = std::hash<double>{}(instance.x);
    std::size_t y_hash = std::hash<double>{}(instance.y);
    return x_hash ^ (y_hash << 1);
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
  Eigen::VectorXd raw_positions;
  std::vector<Vertex> positions;
  // TODO(@jparr721) - Fix xy_hash and xy in general to support 3d.
  std::unordered_map<xy, unsigned int, xy_hash> indices;

  Mesh(const Eigen::VectorXd& vertices, bool is_2d = false);
  Mesh(const Eigen::VectorXd& vertices, const std::vector<QVector3D>& colors);

  void LoadVertices(const Eigen::VectorXd& vertices, bool is_2d);
  void LoadVertices(const Eigen::VectorXd& vertices,
                    const std::vector<QVector3D>& colors);

  void IndexDuplicateVertices();
  void Initialize(QOpenGLBuffer& vbo);
  void Render(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao);

  unsigned int rows() { return raw_positions.rows(); }
  unsigned int index(double x, double y) { return indices.at(xy{x, y}); }
};
