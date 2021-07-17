#pragma once

#include <vector>

#include "BoundaryCondition.h"
#include "Numerics.h"
#include "Rve.h"

#include "Mesh.h"

/**
This class holds the numerical simulator for the corotational linear FEA model.
**/
class LinearTetrahedral {
    using BetaSubMatrixXf = Eigen::Matrix<Real, 6, 3>;
    struct ElementStiffness {
        Matrix12r stiffness_matrix;
        std::vector<int> indices;
    };

  public:
    const unsigned int integrator_size;

    // The boundary forces
    VectorXr boundary_forces;

    // The Per-Element Stiffness Matrix
    MatrixXr per_element_stiffness;

    // The global stiffness matrix
    MatrixXr global_stiffness;

    // The global displacement vector
    VectorXr global_displacement;

    // Element stiffness matrices and mapped coordinates
    std::vector<ElementStiffness> element_stiffnesses;

    // Boundary conditions on nodes in the mesh
    std::vector<BoundaryCondition> boundary_conditions;

    LinearTetrahedral(Real youngs_modulus, Real poissons_ratio,
                      const std::shared_ptr<Mesh>& mesh,
                      std::vector<BoundaryCondition> boundary_conditions);

    void AssembleGlobalStiffness(const std::shared_ptr<Mesh>& mesh);

    /**
    @brief Assemble 12x12 element stiffness matrix. Given by [k] = V[B]^T[D][B]
    where V is the volume of the element
    **/
    void AssembleElementStiffness(Real youngs_modulus, Real poissons_ratio,
                                  const std::shared_ptr<Mesh>& mesh);

    void AssembleBoundaryForces();

    /*
    @brief Calculates the element plane stresses
    */
    MatrixXr AssembleElementPlaneStresses(const MatrixXr& sigmas);

    /*
    @brief Applies the vector of boundary conditions to the nodes and solves the
    dynamic form problem
    */
    MatrixXr Solve(Real youngs_modulus, Real poissons_ratio,
                   const std::shared_ptr<Mesh>& mesh);

    /*
    @brief Applies the vector of boundary conditions to the nodes and solves the
    static form problem via homogenization
    */
    MatrixXr SolveStatic(Real youngs_modulus, Real poissions_ratio,
                         const std::shared_ptr<Mesh>& mesh);

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
    [[nodiscard]] Real ConstructShapeFunctionParameter(Real p1, Real p2,
                                                       Real p3, Real p4,
                                                       Real p5, Real p6);

    [[nodiscard]] VectorXr ComputeRenderedDisplacements(int displacements_size);

    [[nodiscard]] Matrix6r AssembleConstitutiveMatrix(Real youngs_modulus,
                                                      Real poissons_ratio);
    [[nodiscard]] MatrixXr AssembleStrainRelationshipMatrix(
        const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
        const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four);
    [[nodiscard]] MatrixXr
    ComputeElementStress(Real youngs_modulus, Real poissons_ratio,
                         const VectorXr& displacement,
                         const std::shared_ptr<Mesh>& mesh);
};
