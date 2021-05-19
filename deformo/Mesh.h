#pragma once

#include "EigenTypes.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <string>

class Mesh {
 public:
  Eigen::MatrixXf faces;
  Eigen::MatrixXf tetrahedrals;
  Eigen::RowMajorVectorXf positions;
  Eigen::RowMajorVectorXf colors;
  std::vector<unsigned int> indices;

  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXf& F,
       const Eigen::MatrixXf& T);

  void Update();

  /*
  \brief We want the vertex data to be vectorized so that way we can more
  easily do operations on it.

  We also set our static color for rendering.
  */
  void InitializeVertexPositions(const Eigen::MatrixXf& V);

  [[nodiscard]] const float* data() { return positions.data(); }
  [[nodiscard]] int size_bytes() { return positions.size() * sizeof(float); }
  [[nodiscard]] int size() { return positions.size(); }

  [[nodiscard]] const float* colors_data() { return colors.data(); }
  [[nodiscard]] int colors_size_bytes() { return colors.size() * sizeof(float); }
  [[nodiscard]] int colors_size() { return colors.size(); }
};
