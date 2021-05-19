#include "Mesh.h"

#include <Eigen/Core>
#include <iostream>

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXi& F,
           const Eigen::MatrixXi& T) {
  Vectorize(positions, V);
  Vectorize(faces, F);
  Vectorize(tetrahedrals, T);

  colors.resize(positions.rows(), positions.cols());
  for (int i = 0; i < positions.size(); ++i) {
    colors(i) = i % 3 == 0 ? 1.f : 0.f;
  }
};

void Mesh::Vectorize(Eigen::RowMajorVectorXf& out, const Eigen::MatrixXf& in) {
  out.resize(in.rows() * in.cols(), 1);
  Eigen::MatrixXf in_t = in.transpose();
  out = Eigen::Map<Eigen::VectorXf>(in_t.data(), in_t.rows() * in_t.cols());
}

void Mesh::Vectorize(Eigen::RowMajorVectorXi& out, const Eigen::MatrixXi& in) {
  out.resize(in.rows() * in.cols(), 1);
  Eigen::MatrixXi in_t = in.transpose();
  out = Eigen::Map<Eigen::VectorXi>(in_t.data(), in_t.rows() * in_t.cols());
}
