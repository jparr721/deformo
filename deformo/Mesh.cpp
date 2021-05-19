#include "Mesh.h"

#include <Eigen/Core>
#include <iostream>

Mesh::Mesh(const Eigen::MatrixXf& V, const Eigen::MatrixXf& F,
           const Eigen::MatrixXf& T) {
  faces = F.array().colwise() - Eigen::VectorXf::Ones(F.rows()).array();
  tetrahedrals = T.array().colwise() - Eigen::VectorXf::Ones(T.rows()).array();

  InitializeVertexPositions(V);
};

void Mesh::InitializeVertexPositions(const Eigen::MatrixXf& V) {
  //positions = V.transpose();
  positions.resize(V.rows() * V.cols(), 1);
  Eigen::MatrixXf Vt = V.transpose();
  positions = Eigen::Map<Eigen::VectorXf>(Vt.data(), Vt.rows() * Vt.cols());

  indices.reserve(faces.rows());
  for (int i = 0; i < faces.rows(); ++i) {
    const Eigen::Vector3f row = faces.row(i);
    indices.push_back(row.x());
    indices.push_back(row.y());
    indices.push_back(row.z());
  }

  colors.resize(positions.rows(), positions.cols());
  for (int i = 0; i < positions.size(); ++i) {
    colors(i) = i % 3 == 0 ? 1.f : 0.f;
  }
}
