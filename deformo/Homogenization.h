#pragma once

#include "AbstractGenerator.h"
#include "DeformoAssert.h"
#include "Material.h"
#include "Numerics.h"
#include <array>

class Homogenization {
    using MatrixXi = MatrixX<int>;
    using VectorXi = VectorX<int>;

  public:
    Homogenization(const Tensor3r& implicit_surface,
                   const Material& material_1);
    Homogenization(const Tensor3r& implicit_surface, const Material& material_1,
                   const Material& material_2);

    auto E() const noexcept -> Real { return homogenized_E_; }
    auto v() const noexcept -> Real { return homogenized_v_; }
    auto Stiffness() const -> Matrix6r { return constitutive_tensor_; }

    /// <summary>
    /// Solves the integral over the volume of the voxel for the difference of the 
    /// macro and micro scale strain tensors.
    /// </summary>
    auto Solve() -> void;

    /// <summary>
    /// Creates the matrix S by inverting the stiffness constutive tensor C and gets 
    /// the 6x6 compliance matrix which contains our material coefficients.
    /// </summary>
    /// <returns></returns>
    auto ComputeMaterialCoefficients() -> void;

    /// <summary>
    /// Computes the finite element approximation of the stiffness and load
    /// matrices
    /// </summary>
    /// <param name="a"></param>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <returns></returns>
    auto ComputeHexahedron(Real a, Real b, Real c) -> std::array<MatrixXr, 4>;

    auto ComputeElementDegreesOfFreedom(unsigned int n_elements) -> MatrixXi;
    auto ComputeUniqueNodes(unsigned int n_elements) -> Tensor3i;
    auto
    ComputeUniqueDegreesOfFreedom(const MatrixXi& element_degrees_of_freedom,
                                  const Tensor3i& unique_nodes) -> MatrixXi;

    auto AssembleStiffnessMatrix(unsigned int n_degrees_of_freedom,
                                 const MatrixXi& unique_degrees_of_freedom,
                                 const MatrixXr& ke_lambda,
                                 const MatrixXr& ke_mu) -> SparseMatrixXr;
    auto AssembleLoadMatrix(unsigned int n_elements,
                            unsigned int n_degrees_of_freedom,
                            const MatrixXi& unique_degrees_of_freedom,
                            const MatrixXr& fe_lambda, const MatrixXr& fe_mu)
        -> SparseMatrixXr;

    /// <summary>
    /// Solves the finite element method with conjugate gradient with incomplete
    /// cholesky preconditioner.
    /// </summary>
    /// <param name="n_degrees_of_freedom">Total degrees of freedom for all
    /// nodes</param> <param name="stiffness">The stiffness matrix K</param>
    /// <param name="load">The load matrix F</param>
    /// <param name="unique_degrees_of_freedom">The degrees of freedom for
    /// non-void regions</param> 
    /// <returns>Nodal displacement matrix Chi (X_e)</returns>
    auto ComputeDisplacement(unsigned int n_degrees_of_freedom,
                             const MatrixXr& stiffness, const MatrixXr& load,
                             const MatrixXi& unique_degrees_of_freedom)
        -> MatrixXr;
    auto ComputeUnitStrainParameters(unsigned int n_elements,
                                     const std::array<MatrixXr, 4>& hexahedron)
        -> Tensor3r;

  private:
    bool is_one_material_ = false;

    unsigned int cell_len_x_ = 0;
    unsigned int cell_len_y_ = 0;
    unsigned int cell_len_z_ = 0;

    Real homogenized_E_ = 0;
    Real homogenized_v_ = 0;

    Matrix6r constitutive_tensor_;

    Tensor3r lambda_;
    Tensor3r mu_;

    Tensor3r voxel_;

    DeformoAssertion assertion;

    Material primary_material_;

    // Constitutive Tensor Collection
    auto AssembleConstitutiveTensor(const MatrixXi& unique_degrees_of_freedom,
                                    const MatrixXr& ke_lambda,
                                    const MatrixXr& ke_mu,
                                    const MatrixXr& displacement,
                                    const Tensor3r& unit_strain_parameter)
        -> void;
};
