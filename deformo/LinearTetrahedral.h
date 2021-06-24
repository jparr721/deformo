#pragma once

#include <vector>

#include "BoundaryCondition.h"
#include "EigenTypes.h"
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
    const unsigned int integrator_size;

    // The boundary forces
    Eigen::VectorXf boundary_forces;

    // The Per-Element Stiffness Matrix
    Eigen::MatrixXf per_element_stiffness;

    // The global stiffness matrix
    Eigen::MatrixXf global_stiffness;

    // The global displacement vector
    Eigen::VectorXf global_displacement;

    // Element stiffness matrices and mapped coordinates
    std::vector<ElementStiffness> element_stiffnesses;

    // Boundary conditions on nodes in the mesh
    std::vector<BoundaryCondition> boundary_conditions;

    LinearTetrahedral(float youngs_modulus, float poissons_ratio,
                      const std::shared_ptr<Mesh>& mesh,
                      std::vector<BoundaryCondition> boundary_conditions);

    void AssembleGlobalStiffness(const std::shared_ptr<Mesh>& mesh);

    /**
    @brief Assemble 12x12 element stiffness matrix. Given by [k] = V[B]^T[D][B]
    where V is the volume of the element
    **/
    void AssembleElementStiffness(float youngs_modulus, float poissons_ratio, const std::shared_ptr<Mesh>& mesh);

    void AssembleBoundaryForces();

    /*
    @brief Calculates the element plane stresses
    */
    Eigen::MatrixXf AssembleElementPlaneStresses(const Eigen::MatrixXf& sigmas);

    /*
    @brief Applies the vector of boundary conditions to the nodes and solves
    */
    Eigen::MatrixXf Solve(float youngs_modulus, float poissons_ratio, const std::shared_ptr<Mesh>& mesh);

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

    [[nodiscard]] Eigen::VectorXf ComputeRenderedDisplacements(int displacements_size);

    [[nodiscard]] Eigen::Matrix66f
    AssembleStressStrainMatrix(float youngs_modulus,
                               float poissons_ratio);
    [[nodiscard]] Eigen::MatrixXf AssembleStrainRelationshipMatrix(
        const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
        const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four);
};
