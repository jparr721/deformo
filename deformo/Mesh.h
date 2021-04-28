#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

class Mesh {
 public:
  Eigen::VectorXd vertices;
  std::unordered_map<std::array<double, 2>, unsigned int> indices;

  Mesh(const Eigen::VectorXd& vertices_);

  unsigned int rows() { return vertices.rows(); }
  unsigned int cols() { return vertices.cols(); }
  unsigned int index(double x, double y) {
    return indices[std::array<double, 2>{x, y}];
  }

 private:
  unsigned int vao = 0;
  unsigned int vbo = 0;

  void LoadVBO();
};
