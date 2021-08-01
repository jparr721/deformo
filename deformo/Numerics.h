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
template <typename T> using Vector6 = Eigen::Matrix<T, 6, 1>;
template <typename T> using Vector12 = Eigen::Matrix<T, 12, 1>;

template <typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using Matrix2 = Eigen::Matrix<T, 2, 2>;
template <typename T> using Matrix3 = Eigen::Matrix<T, 3, 3>;
template <typename T> using Matrix4 = Eigen::Matrix<T, 4, 4>;

template <typename T>
class Tensor3 {
  public:
    Tensor3() = default;
    Tensor3(const Eigen::Tensor<T, 3> instance) : instance_(instance) {}
    explicit Tensor3(const Vector3<unsigned>& dims) {
        Resize(dims.x(), dims.y(), dims.z());
    }
    Tensor3(unsigned int layers, unsigned int width, unsigned int height) {
        Resize(layers, width, height);
    }

    auto Matrix(unsigned int rows, unsigned int cols) -> MatrixX<T> {
        return Eigen::Map<const MatrixX<T>>(instance_.data(), rows, cols);
    }
    auto Dimension(const unsigned int dim) const -> unsigned int {
        return instance_.dimension(dim);
    }
    auto Dimensions() const -> Vector3<unsigned> {
        return Vector3<unsigned>(Dimension(0), Dimension(1), Dimension(2));
    }
    auto Instance() -> Eigen::Tensor<T, 3>& { return instance_; }

    auto SetConstant(T value) -> void { instance_.setConstant(value); }

    auto Resize(unsigned int layers, unsigned int width, unsigned int height)
        -> void {
        instance_.resize(layers, width, height);
    }

    friend std::ostream& operator<<(std::ostream& out, const Tensor3& t) {
        return out << t.instance_;
    }

    template <typename... Indices>
    void operator()(T value, Indices&&... indices) {
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
    Eigen::Tensor<T, 3> instance_;

    auto Layer(const unsigned int layer) const -> MatrixX<T> {
        MatrixX<T> m;
        const unsigned int width = instance_.dimension(1);
        const unsigned int height = instance_.dimension(2);

        m.resize(width, height);
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                m(i, j) = instance_(i, j, layer);
            }
        }

        return m;
    }

    auto Row(const unsigned int layer, const unsigned int row) const
        -> VectorX<T> {
        VectorX<T> v;
        const unsigned int cols = instance_.dimension(2);
        v.resize(cols);

        for (int col = 0; col < cols; ++col) {
            v(col) = instance_(row, col, layer);
        }

        return v;
    }
};

using Tensor3r = Tensor3<Real>;
using Tensor3i = Tensor3<int>;

namespace linear_algebra {
auto OneDimensionalLinearInterpolation(Real low, Real high, Real interval)
    -> std::vector<Real>;

template <typename T> constexpr auto Lerp(T a, T b, Real t) noexcept -> T {
    return a + (b - a) * t;
};

auto LinSpace(Real start, Real stop, unsigned int num) -> VectorXr;
} // namespace linear_algebra
