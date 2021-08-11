#pragma once

#include "DeformoAssert.h"
#include "Material.h"
#include "Numerics.h"
#include "Rve.h"
#include "Utils.h"
#include <vector>

const auto UnsignedVectorComparison = [](const Vector3<unsigned>& lhs,
                                         const Vector3<unsigned>& rhs) -> bool {
    if (lhs.x() < rhs.x())
        return true;
    if (rhs.x() < lhs.x())
        return false;
    if (lhs.y() < rhs.y())
        return true;
    if (rhs.y() < lhs.y())
        return false;
    return lhs.z() < rhs.z();
};

const auto UnsignedVectorEqual = [](const Vector3<unsigned>& lhs,
                                    const Vector3<unsigned> rhs) -> bool {
    return lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() == rhs.z();
};


template <typename T> class ImplicitSurfaceGenerator {
  public:
    struct Inclusion {
        unsigned int n_inclusions;
        unsigned int area;
        unsigned int depth;
        unsigned int rows;
        unsigned int cols;
    };

    enum class ImplicitSurfaceCharacteristics {
        kIsotropic = 0x00,
        kAnIsotropic = 0x01,
    };

    enum class ImplicitSurfaceMicrostructure {
        kUniform = 0x00,
        kComposite = 0x01,
    };

    int minimum_surface_padding = 1;

    ImplicitSurfaceGenerator(const unsigned int height,
                             const unsigned int width, const unsigned int depth,
                             ImplicitSurfaceCharacteristics behavior,
                             ImplicitSurfaceMicrostructure microstructure,
                             const Inclusion inclusion, const Material& material_1,
                             const Material& material_2)
        : behavior_(behavior), microstructure_(microstructure),
          inclusion_(inclusion), material_1_(material_1),
          material_2_(material_2) {
        implicit_surface_.Resize(height, width, depth);
        implicit_surface_.SetConstant(material_1_.number);
    }

    auto Generate() -> Tensor3<T> {
        return behavior_ == ImplicitSurfaceCharacteristics::kIsotropic
                   ? GenerateIsotropicMaterial()
                   : GenerateAnisotropicMaterial();
    }

  private:
    DeformoAssertion deformo_assert_;

    ImplicitSurfaceCharacteristics behavior_;

    ImplicitSurfaceMicrostructure microstructure_;

    Tensor3<T> implicit_surface_;

    Inclusion inclusion_;

    Material material_1_;
    Material material_2_;

    // Helpers ======================
    auto LayerContainsSecondaryMaterial(const MatrixXr& layer) -> bool {
        const auto rows = layer.rows();
        const auto cols = layer.cols();

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                if (layer(row, col) == material_2_.number) {
                    return true;
                }
            }
        }

        return false;
    }

    /*
    @brief Ensures padding is even on all sides, returning false if not. Only
    applies to isotropic material meshes.
    */
    auto CheckLayerPadding(const VectorXr& layer, int max) -> void {
        const auto left_padding = layer(0) + 1;
        const auto right_padding = max - layer(layer.rows() - 1);

        deformo_assert_.Assert(
            left_padding == right_padding, __FUNCTION__, __FILE__, __LINE__,
            "Padding does not match: ", left_padding, right_padding, layer);
    }

    auto GenerateIsotropicMaterial() -> Tensor3<T> {
        // Y axis origin with padding
        VectorXr y_axis_origins = linear_algebra::LinSpace(
            inclusion_.area + minimum_surface_padding,
            (implicit_surface_.Dimension(1) + 1 -
             (inclusion_.area + minimum_surface_padding)),
            inclusion_.rows);
        y_axis_origins -= VectorXr::Ones(y_axis_origins.rows());
        CheckLayerPadding(y_axis_origins, implicit_surface_.Dimension(1));

        VectorXr x_axis_origins = linear_algebra::LinSpace(
            inclusion_.area + minimum_surface_padding,
            (implicit_surface_.Dimension(0) + 1 -
             (inclusion_.area + minimum_surface_padding)),
            inclusion_.cols);
        x_axis_origins -= VectorXr::Ones(x_axis_origins.rows());
        CheckLayerPadding(x_axis_origins, implicit_surface_.Dimension(0));

        std::vector<Vector2<unsigned int>> centroids;
        for (int i = 0; i < y_axis_origins.rows(); ++i) {
            for (int j = 0; j < x_axis_origins.rows(); ++j) {
                const int y = static_cast<int>(y_axis_origins(i));
                const int x = static_cast<int>(x_axis_origins(j));
                centroids.emplace_back(x, y);
            }
        }

        const VectorXr layer_layout = linear_algebra::LinSpace(
            0, implicit_surface_.Dimension(2) - inclusion_.depth,
            implicit_surface_.Dimension(2) / inclusion_.depth);

        std::vector<Vector3<unsigned int>> indices;
        for (int i = 0; i < layer_layout.rows(); ++i) {
            const int layer = static_cast<int>(layer_layout(i));
            for (const auto& centroid : centroids) {
                const auto cube_indices = MakeShapedIndices(
                    centroid,
                    Vector3<unsigned int>(inclusion_.area, inclusion_.area,
                                          inclusion_.depth),
                    layer);
                indices.insert(indices.end(), cube_indices.begin(),
                               cube_indices.end());
            }
        }

        SetFromIndices(indices);

        return implicit_surface_;
    }

    auto GenerateAnisotropicMaterial() -> Tensor3<T> {
        return implicit_surface_;
    }

    // Shape Generators
    auto MakeShapedIndices(const Vector2<unsigned>& centroid,
                           const Vector3<unsigned>& shape,
                           unsigned int layer_number)
        -> std::vector<Vector3<unsigned int>> {
        std::vector<Vector3<unsigned int>> indices;
        const unsigned int width{shape.x() / 2};
        const unsigned int height{shape.y() / 2};
        const unsigned int depth = shape.z();

        const unsigned int current_x = centroid.x();
        const unsigned int current_y = centroid.y();
        for (auto layer = layer_number; layer < layer_number + depth; ++layer) {
            for (auto x = 0u; x < width; ++x) {
                for (auto y = 0u; y < height; ++y) {
                    indices.emplace_back(current_x + x, current_y + y, layer);
                    indices.emplace_back(current_x - x, current_y - y, layer);
                    indices.emplace_back(current_x - x, current_y + y, layer);
                    indices.emplace_back(current_x + x, current_y - y, layer);
                }
            }
        }

        // TODO(@jparr721) - This is stupid and should be re-written.
        std::sort(indices.begin(), indices.end(), UnsignedVectorComparison);
        indices.erase(
            std::unique(indices.begin(), indices.end(), UnsignedVectorEqual),
            indices.end());

        return indices;
    }

    auto SetFromIndices(const std::vector<Vector3<unsigned int>>& indices)
        -> void {
        for (const Vector3<unsigned int>& index : indices) {
            implicit_surface_(index.x(), index.y(), index.z()) =
                material_2_.number;
        }
    }
};
