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
  positions.resize(V.rows() * V.cols());
  Eigen::MatrixXf Vt = V.transpose();
  positions = Eigen::Map<Eigen::VectorXf>(Vt.data(), Vt.rows() * Vt.cols());

  indices.reserve(faces.rows());
  for (int i = 0; i < faces.rows(); ++i) {
    const Eigen::Vector3f row = faces.row(i);
    indices.push_back(row.x());
    indices.push_back(row.y());
    indices.push_back(row.z());
  }

  // TODO(@jparr721) - Add Permutations for Ele file

  //Eigen::MatrixXf Tt = tetrahedrals.transpose();
  //Eigen::VectorXf tp =
  //    Eigen::Map<Eigen::VectorXf>(Tt.data(), Tt.rows() * Tt.cols());
  //indices.reserve(tp.rows());
  //for (int i = 0; i < tp.rows(); ++i) {
  //  indices.push_back(tp(i));
  //}

  dirty = true;
  Update();
}

void Mesh::Update() {
  // Only update when something changes.
  if (!dirty) {
    return;
  }

  std::vector<Vertex> updated_positions;

  for (int i = 0; i < positions.rows(); i += 3) {
    updated_positions.push_back(
        Vertex{Eigen::Vector3f(positions[i], positions[i + 1], positions[i + 2]),
               Eigen::Vector3f(0.f, 0.f, 1.f)});
  }

  renderable_positions = updated_positions;
  dirty = false;
}