#pragma once

#include <Eigen/Dense>
#include <array>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct xy {
  double x;
  double y;
  bool operator<(const xy& other) const { return x < other.x || (x == other.x && y < other.y); }
  bool operator==(const xy& other) const { return x == other.x && y == other.y; }
};

class Mesh {
 public:
  Eigen::VectorXd vertices;
  std::map<xy, unsigned int> indices;

  Mesh(const Eigen::VectorXd& vertices_);

  unsigned int rows() { return vertices.rows(); }
  unsigned int cols() { return vertices.cols(); }
  unsigned int index(double x, double y) { return indices[xy{x, y}]; }

 private:
  unsigned int vao = 0;
  unsigned int vbo = 0;

  void LoadVBO();
};
