#pragma once

#include <Eigen/Dense>
#include <QVector3D>
#include <array>
#include <functional>
#include <iomanip>
#include <string>
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

  Mesh(const Eigen::VectorXd& vertices) : raw_positions(vertices) {
    LoadVertices(vertices);
    IndexDuplicateVertices();
  }
  Mesh(const Eigen::VectorXd& vertices, const std::vector<QVector3D>& colors)
      : raw_positions(vertices) {
    LoadVertices(vertices, colors);
  }

  void LoadVertices(const Eigen::VectorXd& vertices) {
    for (auto i = 0u; i < raw_positions.rows(); i += 3) {
      positions.push_back(Vertex(QVector3D(
          raw_positions(i), raw_positions(i + 1), raw_positions(i + 2))));
    }
  }

  void LoadVertices(const Eigen::VectorXd& vertices,
                    const std::vector<QVector3D>& colors) {
    for (auto i = 0u, j = 0u; i < raw_positions.rows(); i += 3, ++j) {
      positions.push_back(
          Vertex(QVector3D(raw_positions(i), raw_positions(i + 1),
                           raw_positions(i + 2)),
                 colors[j]));
    }
  }

  void IndexDuplicateVertices() {
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

  void Render() {}

  unsigned int rows() { return raw_positions.rows(); }
  unsigned int index(double x, double y) { return indices.at(xy{x, y}); }
};
