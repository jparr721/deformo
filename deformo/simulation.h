#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EigenTypes.h"

/**
This class holds the numerical simulator for the FEA model.
**/
class Simulation {
 public:
  Simulation(Eigen::Ref<const Eigen::VectorXi> particles_,
             Eigen::Ref<const Eigen::VectorXi> springs_);
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

  // The Mass Matrix
  Eigen::MatrixXd M;

  // The K stiffness matrix
  Eigen::SparseMatrixXd K;

  void AssembleForces();
  void AssembleStiffness();
  void AssembleMassMatrix();
};
