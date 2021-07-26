#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

#ifdef DEFORMO_USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

// Dense Vector Types
using Vector2r = Eigen::Matrix<Real, 2, 1>;
using Vector3r = Eigen::Matrix<Real, 3, 1>;
using Vector4r = Eigen::Matrix<Real, 4, 1>;
using Vector6r = Eigen::Matrix<Real, 6, 1>;
using Vector12r = Eigen::Matrix<Real, 12, 1>;
using VectorXr = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

// Dense Matrix Types
using Matrix2r = Eigen::Matrix<Real, 2, 2>;
using Matrix3r = Eigen::Matrix<Real, 3, 3>;
using Matrix4r = Eigen::Matrix<Real, 4, 4>;
using Matrix6r = Eigen::Matrix<Real, 6, 6>;
using Matrix36r = Eigen::Matrix<Real, 3, 6>;
using Matrix12r = Eigen::Matrix<Real, 12, 12>;
using MatrixXr = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

// Sparse MatrixTypes
using SparseMatrixXr = Eigen::SparseMatrix<Real>;

template <typename T> using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T> using Vector2 = Eigen::Matrix<T, 2, 1>;
template <typename T> using Vector3 = Eigen::Matrix<T, 3, 1>;
template <typename T> using Vector4 = Eigen::Matrix<T, 4, 1>;

template <typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using Matrix2 = Eigen::Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Eigen::Matrix<T, 3, 3>;
template <typename T> using Matrix4 = Eigen::Matrix<T, 4, 4>;

template <unsigned Dim> class Tensor {
  public:
    Tensor() = default;
    Tensor(const Eigen::Tensor<Real, Dim> instance) : instance_(instance) {}
    template <typename... T> explicit Tensor(T&&... dimensions) {
        Resize(dimensions...);
    }

    auto Dimension(const unsigned int dim) const -> unsigned int {
        return instance_.dimension(dim);
    }
    auto Dimensions() const -> unsigned int { return instance_.NumDimensions; }
    auto Instance() -> Eigen::Tensor<Real, Dim>& { return instance_; }

    template <typename... Indices>
    static auto Constant(Real value, Indices&&... indices) -> Tensor<Dim> {
        Tensor<Dim> t;
        t.Resize(indices...);
        constexpr std::size_t n_indices = sizeof...(indices);

        if (n_indices == 3) {
            const auto layers = indices...[0];
            const auto rows = indices...[1];
            const auto cols = indices...[2];

            for (int i = 0; i < layers; ++i) {
                for (int j = 0; j < rows; ++j) {
                    for (int k = 0; k < cols; ++k) {
                        t(value, indices...);
                    }
                }
            }
        }
    }

    template <typename... T> auto Resize(T&&... dimensions) -> void {
        instance_.resize(dimensions...);
    }

    friend std::ostream& operator<<(std::ostream& out, const Tensor& t) {
        return out << t.instance_;
    }

    template <typename... Indices>
    void operator()(Real value, Indices&&... indices) {
        instance_(indices...) = value;
    }

    template <typename... Indices> auto At(Indices&&... indices) const {
        constexpr std::size_t n_indices = sizeof...(indices);

        if constexpr (n_indices == 1) {
            return Layer(indices...);
        }

        if constexpr (n_indices == 2) {
            return Row(indices...);
        }

        if constexpr (n_indices == 3) {
            return instance_(indices...);
        }
    }

  private:
    Eigen::Tensor<Real, Dim> instance_;

    auto Layer(const unsigned int layer) const -> MatrixXr {
        MatrixXr m;
        const unsigned int width = instance_.dimension(1);
        const unsigned int height = instance_.dimension(2);

        m.resize(width, height);
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                m(i, j) = instance_(layer, i, j);
            }
        }

        return m;
    }

    auto Row(const unsigned int layer, const unsigned int row) const
        -> VectorXr {
        VectorXr v;
        const unsigned int cols = instance_.dimension(2);
        v.resize(cols);

        for (int col = 0; col < cols; ++col) {
            v(col) = instance_(layer, row, col);
        }

        return v;
    }
};

// Dense Tensor Types
using Tensor3r = Tensor<3>;

namespace linear_algebra {
auto OneDimensionalLinearInterpolation(Real low, Real high, Real interval)
    -> std::vector<Real>;

template <typename T> constexpr auto Lerp(T a, T b, Real t) noexcept -> T {
    return a + (b - a) * t;
};

auto LinSpace(Real start, Real stop, unsigned int num) -> VectorXr;
} // namespace linear_algebra
