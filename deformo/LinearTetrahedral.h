#pragma once

#include <Eigen/Sparse>
#include <vector>

#include "BoundaryCondition.h"
#include "EigenTypes.h"
#include "Integrators.h"
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
    constexpr static int kStride = 12;
    constexpr static int kFaceStride = 4;

    // \brief Thickness
    constexpr static float kThickness = 2.5e-2f;

    // Timestep constants
    float current_time = 0.f;
    float timestep_size = 0.005f;

    // Modulus of Elasticity
    const float kModulusOfElasticity;

    // Poisson's Ratio
    const float kPoissonsRatio;

    // Integration constants
    float a0 = 1e-10f;
    float a1 = 1e-10f;
    float a2 = 1e-10f;
    float a3 = 1e-10f;

    // The Global Force Vector (Interaction Forces)
    Eigen::VectorXf global_force;

    // Global stacked vertex vector (xy) with index mapping of node orientations
    std::shared_ptr<Mesh> mesh;

    // Global stacked acceleration vector
    Eigen::VectorXf acceleration;

    // Global stacked velocity vector
    Eigen::VectorXf velocity;

    // The Mass Matrix
    Eigen::SparseMatrixXf mass;

    // The Effective Mass Matrix
    Eigen::FullPivLU<Eigen::MatrixXf> effective_mass;

    // The global stiffness matrix
    Eigen::MatrixXf global_stiffness;

    // Element stiffness matrices and mapped coordinates
    std::vector<ElementStiffness> element_stiffnesses;

    // Element Stress vectors for each group of points
    std::vector<Eigen::Vector6f> sigmas;

    // The Element Stresses from each stress vector
    std::vector<Eigen::Vector3f> element_stresses;

    // Boundary conditions on nodes in the mesh
    std::vector<BoundaryCondition> boundary_conditions;

    // The integrator for our sim
    std::unique_ptr<ExplicitCentralDifferenceMethod> integrator;

    LinearTetrahedral(const float modulus_of_elasticity,
                      const float poissons_ratio, const float point_mass,
                      std::shared_ptr<Mesh> mesh,
                      std::vector<BoundaryCondition> boundary_conditions);
    void Update();
    void Integrate();

    void AssembleForces();
    void AssembleMassMatrix(const float point_mass);

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
    void AssembleElementStresses(const Eigen::VectorXf& u,
                                 const Eigen::MatrixXf& B);

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
    float ConstructShapeFunctionParameter(float p1, float p2, float p3,
                                          float p4, float p5, float p6);

    void InitializeIntegrator();

  private:
    // constexpr static unsigned int stride = 3 * kTetrahedronElementCount;

    void InitializeVelocity();
    void InitializeAcceleration();
    void InitializeIntegrationConstants();
    Eigen::VectorXf SolveU(const Eigen::MatrixXf& k, const Eigen::VectorXf& f,
                           const Eigen::VectorXi& indices);
};
