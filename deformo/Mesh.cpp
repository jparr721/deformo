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
  indices = std::vector<float>(vertices.rows());
  std::iota(indices.begin(), indices.end(), 0);

  dirty = true;
  Update();
}

void Mesh::Update() {
  // Only update when something changes.
  if (!dirty) {
    return;
  }

  std::vector<Vertex> updated_positions;

  for (int i = 0; i < vertices.rows(); i += 3) {
    updated_positions.push_back(
        Vertex{QVector3D(vertices[i], vertices[i + 1], vertices[i + 2]),
               QVector3D(0.f, 0.f, 1.f)});
  }

  positions = updated_positions;
  dirty = false;
}