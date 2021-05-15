#pragma once

#include <Eigen/Dense>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <memory>
#include <string>

class Mesh {
 public:
  Eigen::VectorXf colors;
  Eigen::VectorXf vertices;
  const Eigen::MatrixXf faces;
  const Eigen::MatrixXf tetrahedrals;

  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXf& F,
       const Eigen::MatrixXf& T);

  /*
  \brief We want the vertex data to be vectorized so that way we can more
  easily do operations on it.

  We also set our static color for rendering.
  */
  void InitializeVertexPositions(const Eigen::MatrixXf& V);

  [[nodiscard]] int size_bytes() { return vertices.size() * sizeof(float); }
  [[nodiscard]] int size() { return vertices.size(); }
};
