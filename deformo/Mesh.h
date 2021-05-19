#pragma once

#include "EigenTypes.h"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <string>

class Mesh {
 public:
  Eigen::RowMajorVectorXi faces;
  Eigen::RowMajorVectorXi tetrahedrals;
  Eigen::RowMajorVectorXf positions;
  Eigen::RowMajorVectorXf colors;

  Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
       const Eigen::MatrixXi& T);

  [[nodiscard]] const float* data() { return positions.data(); }
  [[nodiscard]] int size_bytes() { return positions.size() * sizeof(float); }
  [[nodiscard]] int size() { return positions.size(); }

  [[nodiscard]] const float* colors_data() { return colors.data(); }
  [[nodiscard]] int colors_size_bytes() { return colors.size() * sizeof(float); }
  [[nodiscard]] int colors_size() { return colors.size(); }

  [[nodiscard]] const int* faces_data() { return faces.data(); }
  [[nodiscard]] int faces_size_bytes() { return faces.size() * sizeof(float); }
  [[nodiscard]] int faces_size() { return faces.size(); }

  private:
  void Vectorize(Eigen::RowMajorVectorXf& out, const Eigen::MatrixXf& in);
  void Vectorize(Eigen::RowMajorVectorXi& out, const Eigen::MatrixXi& in);
};
