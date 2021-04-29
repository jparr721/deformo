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

class Mesh {
 public:
  Eigen::VectorXd vertices;
  std::unordered_map<xy, unsigned int, xy_hash> indices;

  Mesh(const Eigen::VectorXd& vertices_);

  unsigned int rows() { return vertices.rows(); }
  unsigned int cols() { return vertices.cols(); }
  unsigned int index(double x, double y) { return indices.at(xy{x, y}); }

 private:
  unsigned int vao = 0;
  unsigned int vbo = 0;

  void LoadVBO();
};
