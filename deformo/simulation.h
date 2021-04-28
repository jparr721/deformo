#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <utility>
#include <vector>

#include "EigenTypes.h"
#include "Mesh.h"

/**
This class holds the numerical simulator for the corotational linear FEA model.
**/
class Simulation {
 private:
  struct ElementStiffness {
    Eigen::Matrix66d stiffness_matrix;
    std::vector<unsigned int> indices;
  };

 public:
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

  // Integration constants
  double a0 = 1e-10;
  double a1 = 1e-10;
  double a2 = 1e-10;
  double a3 = 1e-10;

  // The Global Force Vector (Interaction Forces)
  Eigen::VectorXd F_ext;

  // The Global Load Vector (at current_time)
  Eigen::VectorXd R_hat;

  // Global stacked vertex vector from the last timestep
  Eigen::VectorXd last_displacement;

  // Global stacked vertex vector (xy)
  std::shared_ptr<Mesh> mesh;

  // Global stacked acceleration vector
  Eigen::VectorXd acceleration;

  // Global stacked velocity vector
  Eigen::VectorXd velocity;

  // The Mass Matrix
  Eigen::SparseMatrixXd M;

  // The Effective Mass Matrix
  Eigen::SparseMatrixXd M_hat;

  // The solved version of M_hat
  Eigen::SparseMatrixXd M_hat_inverse;

  // The global stiffness matrix
  Eigen::MatrixXd K;

  // Element stiffness matrices and mapped coordinates
  std::vector<ElementStiffness> k;

  Simulation() = default;
  Simulation(double mass_, double E_, double NU_,
             const std::shared_ptr<Mesh>& mesh_);
  void Update();
  void Integrate();

  void AssembleForces();
  void AssembleGlobalStiffness();
  void AssembleElementStiffness();
  void AssembleMassMatrix();

  void InitializeVelocity();
  void InitializeAcceleration();
  void InitializeIntegrationConstants();
};
