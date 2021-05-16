#pragma once

#include <Eigen/Dense>
#include <QVector3D>
#include <vector>
#include <memory>
#include <string>

struct Vertex {
  QVector3D position;
  QVector3D color;
};

class Mesh {
 public:
  bool dirty = false;

  Eigen::MatrixXf faces;
  Eigen::MatrixXf tetrahedrals;
  Eigen::VectorXf positions;
  std::vector<Vertex> renderable_positions;
  std::vector<unsigned short> indices;

  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXf& F,
       const Eigen::MatrixXf& T);

  void Update();

  /*
  \brief We want the vertex data to be vectorized so that way we can more
  easily do operations on it.

  We also set our static color for rendering.
  */
  void InitializeVertexPositions(const Eigen::MatrixXf& V);

  [[nodiscard]] const Vertex* data() { return renderable_positions.data(); }
  [[nodiscard]] int size_bytes() { return renderable_positions.size() * sizeof(Vertex); }
  [[nodiscard]] int size() { return renderable_positions.size(); }
  //[[nodiscard]] int size_bytes() { return vertices.size() * sizeof(float); }
  //[[nodiscard]] int size() { return vertices.size(); }
};
