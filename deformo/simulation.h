#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "EigenTypes.h"

/**
This class holds the numerical simulator for the FEA model. This system assumes
uniform springs density based on the particle positions.
**/
class Simulation {
 public:
  Simulation(Eigen::Ref<const Eigen::VectorXd> particles_, const std::vector<Eigen::Vector3d>& fixed_indices);
  void Update();
  void Integrate();

 private:
  // Our starting rest length of the two particles
  int rest_length = 0;

  // Timestep constants
  double current_time = 0.;
  double timestep_size = 0.005;

  // Spring stiffness constant
  double k = 100;

  // Point mass
  double point_mass = 1.;

  // The Global Force Vector
  Eigen::VectorXd F_ext;

  // Global stacked particle vector (xyz)
  Eigen::VectorXd particles;

  // Global stacked spring vector (xyz)
  Eigen::VectorXd springs;

  // Selection Matrix for particles
  Eigen::SparseMatrixXd P;

  // The Mass Matrix
  Eigen::MatrixXd M;

  // The K stiffness matrix
  Eigen::SparseMatrixXd K;

  void AssembleForces();
  void AssembleStiffness();
  void AssembleMassMatrix();
  void AssembleFixedPointSelectionMatrices(const std::vector<Eigen::Vector3d>& fixed_indices);
};
