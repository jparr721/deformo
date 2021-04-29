#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <utility>
#include <vector>

#include "EigenTypes.h"
#include "Mesh.h"

/*
Static boundary condition for a given node
*/
struct BoundaryCondition {
  unsigned int node;
  xy force;
};

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

  // Global stacked vertex vector (xy) with index mapping of node orientations
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

  // The global displacement vector
  Eigen::VectorXd U;

  // The nodal displacement vector
  std::vector<Eigen::Vector6d> nodal_displacements;

  // Element stiffness matrices and mapped coordinates
  std::vector<ElementStiffness> k;

  // Element Stress vectors for each group of points
  std::vector<Eigen::Vector3d> sigmas;

  // The Plane Stresses from each stress vector
  std::vector<Eigen::Vector3d> plane_stresses;

  // Boundary conditions on nodes in the mesh
  std::vector<BoundaryCondition> boundary_conditions;

  Simulation() = default;
  Simulation(double mass_, double E_, double NU_,
             const std::shared_ptr<Mesh>& mesh_,
             const std::vector<BoundaryCondition>& boundary_conditions_);
  void Update();
  void Integrate();

  void AssembleForces();
  void AssembleMassMatrix();

  void AssembleGlobalStiffness();

  void AssembleStressStrainMatrix();
  void AssembleStrainRelationshipMatrix(double xi, double xj, double xm,
                                        double yi, double yj, double ym);

  /**
  @brief Calculates the per-element stiffness matrix
  **/
  void AssembleElementStiffness();

  /*
  @brief Calculates the per-element stresses using our tensile parameters
  */
  void AssembleElementStresses();

  /*
  @bried Calculates the element principal stresses
  */
  void AssembleElementPlaneStresses();

  /*
  @brief Applies the vector of boundary conditions to the nodes and solves
  */
  void Solve();

 private:
  // Strain Relationship Matrix
  Eigen::Matrix36d B;

  // Stress-Strain Matrix
  Eigen::Matrix3d D;

  void InitializeVelocity();
  void InitializeAcceleration();
  void InitializeIntegrationConstants();
  void SolveU(Eigen::MatrixXd k, Eigen::VectorXd f);
};
