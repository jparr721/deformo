#include "Mesh.h"

#include <Eigen/Core>
#include <iostream>

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
           const Eigen::MatrixXi& T) {
  Vectorize<Eigen::VectorXf>(positions, V);
  Vectorize<Eigen::VectorXi>(faces, F);
  Vectorize<Eigen::VectorXi>(tetrahedrals, T);

  colors.resize(positions.rows(), positions.cols());
  for (int i = 0; i < positions.size(); ++i) {
    colors(i) = i % 3 == 0 ? 1.f : 0.f;
  }
};
