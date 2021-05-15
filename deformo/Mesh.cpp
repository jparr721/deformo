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
  vertices.resize(V.rows() * V.cols());
  Eigen::MatrixXf Vt = V.transpose();
  vertices = Eigen::Map<Eigen::VectorXf>(Vt.data(), Vt.rows() * Vt.cols());

  colors.resize(vertices.rows());
  for (int i = 0; i < vertices.rows(); i += 3) {
    colors.segment(i, 3) << 1.0f, 0.f, 0.f;
    positions.push_back(
        Vertex{QVector3D(vertices[i], vertices[i + 1], vertices[i + 2]),
               QVector3D(0.f, 0.f, 1.f)});
  }
}