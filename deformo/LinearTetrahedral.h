#pragma once

#include <Eigen/Sparse>
#include <vector>

#include "BoundaryCondition.h"
#include "EigenTypes.h"
#include "ExplicitCentralDifference.h"
#include "Mesh.h"

/**
This class holds the numerical simulator for the corotational linear FEA model.
**/
class LinearTetrahedral {
    using BetaSubMatrixXf = Eigen::Matrix<float, 6, 3>;
    struct ElementStiffness {
        Eigen::Matrix12f stiffness_matrix;
        std::vector<int> indices;
    };

  public:
    constexpr static int kFaceStride = 4;

    const unsigned int integrator_size;

    // Timestep constants
    float current_time = 0.f;
    float dt = 0.1f;

    // Modulus of Elasticity
    const float modulus_of_elasticity;

    // Poisson's Ratio
    const float poissons_ratio;

    // The boundary force indices
    Eigen::VectorXi boundary_force_indices;

    // The boundary forces
    Eigen::VectorXf boundary_forces;

    // Global stacked vertex vector (xy) with index mapping of node orientations
    std::shared_ptr<Mesh> mesh;

    // The Mass Matrix
    Eigen::SparseMatrixXf mass;

    // The Per-Element Stiffness Matrix
    Eigen::MatrixXf per_element_stiffness;

    // The global stiffness matrix
    Eigen::MatrixXf global_stiffness;

    // The global displacement vector
    Eigen::VectorXf global_displacement;

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

    void InitializeIntegrator();

    void AssembleMassMatrix(const float point_mass);

    void AssembleGlobalStiffness();

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

    void AssembleBoundaryForces();

    /*
    @bried Calculates the element plane stresses
    */
    void AssembleElementPlaneStresses();

    /*
    @brief Applies the vector of boundary conditions to the nodes and solves
    */
    void Solve();

    /*
    @brief Compute the volume of the tetrahedral element.
    */
    [[nodiscard]] float ComputeElementVolume(const Eigen::Vector3f& shape_one,
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
    [[nodiscard]] float ConstructShapeFunctionParameter(float p1, float p2,
                                                        float p3, float p4,
                                                        float p5, float p6);

    [[nodiscard]] Eigen::VectorXf ComputeRenderedDisplacements();

    [[nodiscard]] Eigen::Matrix66f
    AssembleStressStrainMatrix(float poissons_ratio,
                               float modulus_of_elasticity);
    [[nodiscard]] Eigen::MatrixXf AssembleStrainRelationshipMatrix(
        const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
        const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four);
};
