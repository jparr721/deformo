#pragma once

#include "DeformoAssert.h"
#include "Numerics.h"
#include "Rve.h"
#include <array>

class Homogenization {
    using MatrixXi = MatrixX<unsigned int>;

  public:
    explicit Homogenization(std::shared_ptr<Rve> rve);

    auto ComputeHexahedron(Real a, Real b, Real c) -> std::array<MatrixXr, 4>;

    auto ComputeDegreesOfFreedom(unsigned int n_elements) -> MatrixX<int>;
    auto ComputeUniqueNodes(unsigned int n_elements) -> MatrixX<int>;

    // Stiffness Calculations
    auto AssembleStiffnessMatrix(const MatrixXi& element_degrees_of_freedom,
                                 const MatrixXr& stiffness_lambda,
                                 const MatrixXr& stiffness_mu,
                                 unsigned int n_degrees_of_freedom) -> MatrixXr;
    auto AssembleLoadMatrix(const MatrixXi& element_degrees_of_freedom,
                            const MatrixXr& stiffness_lambda,
                            const MatrixXr& stiffness_mu,
                            unsigned int n_degrees_of_freedom) -> MatrixXr;
    auto
    ComputeUniqueDegreesOfFreedom(const MatrixXi& element_degrees_of_freedom,
                                  const MatrixXi& unique_nodes) -> MatrixXi;
    auto ComputeDisplacement(const MatrixXr& stiffness, const MatrixXr& load,
                             const MatrixXi& element_degrees_of_freedom,
                             unsigned int n_degrees_of_freedom) -> MatrixXr;
    auto ComputeUnitStrainParameters() -> MatrixXr;

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

    // Utilities
};
