#pragma once

#include "DeformoAssert.h"
#include "Numerics.h"
#include "Rve.h"
#include <array>

class Homogenization {
    using MatrixXi = MatrixX<int>;
    using VectorXi = VectorX<int>;

  public:
    explicit Homogenization(std::shared_ptr<Rve> rve);

    auto E() -> Real { return homogenized_E_; }
    auto v() -> Real { return homogenized_v_; }

    auto Solve() -> void;

    auto ComputeHexahedron(Real a, Real b, Real c) -> std::array<MatrixXr, 4>;

    auto ComputeElementDegreesOfFreedom(unsigned int n_elements) -> MatrixXi;
    auto ComputeUniqueNodes(unsigned int n_elements) -> Tensor3i;
    auto
    ComputeUniqueDegreesOfFreedom(const MatrixXi& element_degrees_of_freedom,
                                  const Tensor3i& unique_nodes) -> MatrixXi;

    // Stiffness Calculations
    auto AssembleStiffnessMatrix(unsigned int n_degrees_of_freedom,
                                 const MatrixXi& unique_degrees_of_freedom,
                                 const MatrixXr& ke_lambda,
                                 const MatrixXr& ke_mu) -> SparseMatrixXr;
    auto AssembleLoadMatrix(unsigned int n_elements,
                            unsigned int n_degrees_of_freedom,
                            const MatrixXi& unique_degrees_of_freedom,
                            const MatrixXr& fe_lambda, const MatrixXr& fe_mu)
        -> SparseMatrixXr;
    auto ComputeDisplacement(unsigned int n_degrees_of_freedom,
                             const MatrixXr& stiffness, const MatrixXr& load,
                             const MatrixXi& element_degrees_of_freedom)
        -> MatrixXr;
    auto ComputeUnitStrainParameters(unsigned int n_elements,
                                     const std::array<MatrixXr, 4>& hexahedron)
        -> Tensor3r;

  private:
    unsigned int cell_len_x_;
    unsigned int cell_len_y_;
    unsigned int cell_len_z_;

    Real homogenized_E_;
    Real homogenized_v_;

    Matrix6r constitutive_tensor_;

    Tensor3r lambda_;
    Tensor3r mu_;

    Tensor3r voxel_;

    std::shared_ptr<Rve> rve_;

    DeformoAssertion assertion;
};
