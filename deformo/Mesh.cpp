#include "Mesh.h"

#include <Eigen/Core>

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXf& F,
           const Eigen::MatrixXf& T)
    : faces(F), tetrahedrals(T) {
  InitializeVertexPositions(V);
};

void Mesh::InitializeVertexPositions(const Eigen::MatrixXf& V) {
  vertices.resize(V.rows() * V.cols());
  Eigen::MatrixXf Vt = V.transpose();
  vertices = Eigen::Map<Eigen::VectorXf>(Vt.data(), Vt.rows() * Vt.cols());

  colors.resize(vertices.rows());
  for (int i = 0; i < vertices.rows(); i += 3) {
    colors.segment(i, 3) << 1.0f, 0.f, 0.f;
  }
}
