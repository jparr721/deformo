#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <vector>
#include <utility>

#include "EigenTypes.h"

/**
This class holds the numerical simulator for the corotational linear FEA model.
**/
class Simulation {
 public:
  Simulation(double mass_, double damping_, double E_, double NU_,
             Eigen::Ref<const Eigen::VectorXd> displacements_);
  void Update();
  void Integrate();

 private:
  constexpr static int PLANE_STRESS = 1;
  constexpr static int PLANE_STRAIN = 2;

  // Timestep constants
  double current_time = 0.;
  double timestep_size = 0.005;

  // Modulus of Elasticity
  double E;

  // Poisson's Ratio
  double NU;

  // Thickness
  double t = 2.5e-2;

  // Point mass
  double mass = 1.;

  // Damping
  double damping = 1.;

  // Integration constants
  double a0 = 1e-10;
  double a1 = 1e-10;

  // The Global Force Vector
  Eigen::VectorXd F_ext;

  // Global stacked vertex vector from the last timestep
  Eigen::VectorXd last_displacement;

  // Global stacked vertex vector (xy)
  Eigen::VectorXd displacements;

  // Global stacked acceleration vector
  Eigen::VectorXd acceleration;

  // Global stacked velocity vector
  Eigen::VectorXd velocity;

  // The Mass Matrix
  Eigen::SparseMatrixXd M;

  // The Effective Mass Matrix
  Eigen::SparseMatrixXd M_hat;

  // The global stiffness matrix
  Eigen::MatrixXd K;

  // Element stiffness matrices and their node numbers
  std::vector<std::pair<Eigen::Matrix66d, std::vector<uint64_t>>> k;

  void AssembleForces();
  void AssembleGlobalStiffness();
  void AssembleElementStiffness();
  void AssembleMassMatrix();

  void InitializeVelocity();
  void InitializeAcceleration();
  void InitializeIntegrationConstants();
};
