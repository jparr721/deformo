#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BoundaryCondition.h"
#include "EigenTypes.h"
#include "Mesh.h"

/**
This class holds the numerical simulator for the corotational linear FEA model.
**/
class LinearTetrahedral {
 private:
  using BetaSubmatrixf = Eigen::Matrix<float, 6, 3>;
  struct ElementStiffness {
    Eigen::Matrix12f stiffness_matrix;
    std::vector<int> indices;
  };

 public:
  constexpr static int PLANE_STRESS = 1;
  constexpr static int PLANE_STRAIN = 2;
  constexpr static int kStride = 12;
  constexpr static int kFaceStride = 4;

  // Timestep constants
  float current_time = 0.;
  float timestep_size = 0.005;

  // Modulus of Elasticity
  float E;

  // Poisson's Ratio
  float NU;

  // Thickness
  float t = 2.5e-2;

  // Point mass
  float mass = 1.;

  // Integration constants
  float a0 = 1e-10;
  float a1 = 1e-10;
  float a2 = 1e-10;
  float a3 = 1e-10;

  // The Global Force Vector (Interaction Forces)
  Eigen::VectorXf F_ext;

  // The Global Load Vector (at current_time)
  Eigen::VectorXf R_hat;

  // Global stacked vertex vector from the last timestep
  Eigen::VectorXf last_displacement;

  // Global stacked vertex vector (xy) with index mapping of node orientations
  std::shared_ptr<Mesh> mesh;

  // Global stacked acceleration vector
  Eigen::VectorXf acceleration;

  // Global stacked velocity vector
  Eigen::VectorXf velocity;

  // The Mass Matrix
  Eigen::SparseMatrixXf M;

  // The Effective Mass Matrix
  Eigen::SparseMatrixXf M_hat;

  // The global stiffness matrix
  Eigen::MatrixXf K;

  // Element stiffness matrices and mapped coordinates
  std::vector<ElementStiffness> k;

  // Element Stress vectors for each group of points
  std::vector<Eigen::Vector6f> sigmas;

  // The Plane Stresses from each stress vector
  std::vector<Eigen::Vector3f> plane_stresses;

  // Boundary conditions on nodes in the mesh
  std::vector<BoundaryCondition> boundary_conditions;

  LinearTetrahedral() = default;
  LinearTetrahedral(float mass_, float E_, float NU_,
                    std::shared_ptr<Mesh>& mesh_,
                    const std::vector<BoundaryCondition>& boundary_conditions_);
  void Update();
  void Integrate();

  void AssembleForces();
  void AssembleMassMatrix();

  void AssembleGlobalStiffness();

  void AssembleStressStrainMatrix(Eigen::Matrix66f& D);
  void AssembleStrainRelationshipMatrix(Eigen::MatrixXf& strain_relationship,
                                        const Eigen::Vector3f& shape_one,
                                        const Eigen::Vector3f& shape_two,
                                        const Eigen::Vector3f& shape_three,
                                        const Eigen::Vector3f& shape_four);

  /**
  @brief Assemble 12x12 element stiffness matrix. Given by [k] = V[B]^T[D][B]
  where V is the volume of the element
  **/
  void AssembleElementStiffness();

  void ComputeElementStiffness(Eigen::Matrix12f& element_stiffness,
                               const Eigen::Vector3f& shape_one,
                               const Eigen::Vector3f& shape_two,
                               const Eigen::Vector3f& shape_three,
                               const Eigen::Vector3f& shape_four);

  /*
  @brief Calculates the per-element stresses using our tensile parameters
  */
  void AssembleElementStresses(Eigen::VectorXf nodal_displacement);

  /*
  @bried Calculates the element principal stresses
  */
  void AssembleElementPlaneStresses();

  /*
  @brief Applies the vector of boundary conditions to the nodes and solves
  */
  void Solve();

  /*
  @brief Compute the volume of the tetrahedral element.
  */
  float ComputeElementVolume(const Eigen::Vector3f& shape_one,
                             const Eigen::Vector3f& shape_two,
                             const Eigen::Vector3f& shape_three,
                             const Eigen::Vector3f& shape_four);

  /*
  @brief Construct the shape function parameter matrix determinant.
  The parameters have pretty terrible naming, they are as follows
  @param p1 Top row, value 1
  @param p2 Top row, value 2
  @param p3 Mid row, value 1
  @param p4 Mid row, value 2
  @param p5 Bot row, value 1
  @param p6 Bot row, value 2
  */
  float ConstructShapeFunctionParameter(float p1, float p2, float p3, float p4,
                                        float p5, float p6);

 private:
  // constexpr static unsigned int stride = 3 * kTetrahedronElementCount;

  void InitializeVelocity();
  void InitializeAcceleration();
  void InitializeIntegrationConstants();
  Eigen::VectorXf SolveU(Eigen::MatrixXf k, Eigen::VectorXf f,
                         Eigen::VectorXi indices);
};
