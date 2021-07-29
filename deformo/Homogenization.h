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
    /*
    \brief Flatten in fortran order.
    */
    template <typename T>
    auto Flatten(const MatrixX<T>& value) -> VectorX<T> {
        // TODO(@jparr721) - Can probably use map
        const auto shape = value.rows() * value.cols();
        return VectorX<T>(Eigen::Map<VectorX<T>> value.data(), shape);
    }

    auto Where(const Tensor3r& input, unsigned int value) const -> Tensor3r {
        const auto layers = input.Dimension(0);
        const auto rows = input.Dimension(1);
        const auto cols = input.Dimension(2);
        Tensor3r output(layers, rows, cols);
        output.SetConstant(1);

        for (auto i = 0u; i < layers; ++i) {
            for (auto j = 0u; j < rows; ++j) {
                for (auto k = 0u; k < cols; ++k) {
                    if (input.At(i, j, k) != value) {
                        output(0, i, j, k);
                    }
                }
            }
        }

        return output;
    }
};
