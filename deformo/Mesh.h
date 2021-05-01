#pragma once

#include <Eigen/Dense>
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

template <unsigned dim>
struct Vertex {
  using V = Eigen::Matrix<double, dim, 1>;

  V position;
  Eigen::Vector3d color;

  Vertex(const V& vertex_)
      : position(vertex_), color(Eigen::Vector3d(1., 0., 0.)) {}
  Vertex(const V& vertex_, const Eigen::Vector3d& color_)
      : position(vertex_), color(color_) {}

  static inline int PositionOffset() { return offsetof(Vertex, position); }
  static inline int ColorOffset() { return offsetof(Vertex, color); }
  static inline int Size() { return dim; }
  static inline int Stride() { return sizeof(position); }
};

template <unsigned dim>
class Mesh {
  using V = Eigen::Matrix<double, dim, 1>;

 public:
  Eigen::VectorXd raw_positions;
  std::vector<Vertex<dim>> positions;
  std::unordered_map<xy, unsigned int, xy_hash> indices;

  Mesh(const Eigen::VectorXd& vertices) : raw_positions(vertices) {
    LoadVertices(vertices);
    IndexDuplicateVertices();
  }
  Mesh(const Eigen::VectorXd& vertices,
       const std::vector<Eigen::Vector3d>& colors)
      : raw_positions(vertices) {
    LoadVertices(vertices, colors);
  }

  void LoadVertices(const Eigen::VectorXd& vertices) {
    for (auto i = 0u; i < raw_positions.rows(); i += 2) {
      positions.push_back(Vertex(V(raw_positions(i), raw_positions(i + 2))));
    }
  }

  void LoadVertices(const Eigen::VectorXd& vertices,
                    const std::vector<Eigen::Vector3d>& colors) {
    int j = 0;
    for (auto i = 0u; i < raw_positions.rows(); i += 2) {
      positions.push_back(Vertex(V(raw_positions(i), raw_positions(i + 2)), colors(j));
		++j;
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

  void Render() {

  }

  unsigned int rows() { return raw_positions.rows(); }
  unsigned int index(double x, double y) { return indices.at(xy{x, y}); }
};
