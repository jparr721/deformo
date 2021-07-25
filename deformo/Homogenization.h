#pragma once

#include "Numerics.h"
#include "Rve.h"

class Homogenization {
    using MatrixXi = MatrixX<unsigned int>;

  public:
    Homogenization(unsigned int cell_len_x, unsigned int cell_len_y,
                   unsigned int cell_len_z, Vector2r lambda, Vector2r mu,
                   std::shared_ptr<Rve> rve);

  private:
    unsigned int cell_len_x_;
    unsigned int cell_len_y_;
    unsigned int cell_len_z_;

    Real lambda_;
    Real mu_;

    Real homogenized_E_;
    Real homogenized_v_;

    Matrix6r constitutive_tensor;

    std::shared_ptr<Rve> rve_;

    auto ComputeDegreesOfFreedom(unsigned int n_elements) -> MatrixXr;
    auto ComputeUniqueNodes(unsigned int n_elements) -> MatrixXr;

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
    auto ComputeDisplacement(const MatrixXr& stiffness,
                             const MatrixXr& load,
                             const MatrixXi& element_degrees_of_freedom,
                             unsigned int n_degrees_of_freedom) -> MatrixXr;
    auto ComputeUnitStrainParameters() -> MatrixXr;

    // Utilities
    template <typename Derived>
    auto Flat2d(Eigen::PlainObjectBase<Derived>& value) -> void;
    template <typename Derived>
    auto Flat1d(Eigen::PlainObjectBase<Derived>& value) -> void;
};
