#include "ImplicitSurfaceGenerator.h"
#include "Utils.h"

//
namespace {
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
} // namespace

ImplicitSurfaceGenerator::ImplicitSurfaceGenerator(const unsigned int height,
                                                   const unsigned int width,
                                                   const unsigned int depth) {
    implicit_surface.Resize(height, width, depth);
    implicit_surface.Instance().setConstant(BinaryMaterial::kPrimaryMaterial);
}

Tensor3r ImplicitSurfaceGenerator::Generate(const BinaryInclusion inclusion) {
    // Y axis origin with padding
    const Vector3r y_axis_origins =
        linear_algebra::LinSpace(inclusion.area + minimum_surface_padding,
                                 implicit_surface.Dimension(2) -
                                     (inclusion.area + minimum_surface_padding),
                                 inclusion.rows);

    std::vector<Vector2<unsigned int>> centroids;
    for (int row = 0; row < inclusion.rows; ++row) {
        const VectorXr x_axis_layout = linear_algebra::LinSpace(
            inclusion.area + minimum_surface_padding,
            implicit_surface.Dimension(2) -
                (inclusion.area + minimum_surface_padding),
            inclusion.cols);

        // Yes, this _is_ prone to rounding error...
        const int y_axis_value = static_cast<int>(y_axis_origins(row));

        for (int i = 0; i < x_axis_layout.rows(); ++i) {
            centroids.emplace_back(x_axis_layout(i), y_axis_value);
        }
    }

    const VectorXr layer_layout = linear_algebra::LinSpace(
        0, implicit_surface.Dimension(0) - inclusion.depth,
        implicit_surface.Dimension(0) / inclusion.depth);

    std::vector<Vector3<unsigned int>> indices;
    for (int i = 0; i < layer_layout.rows(); ++i) {
        const int layer = static_cast<int>(layer_layout(i));
        for (const auto& centroid : centroids) {
            const auto cube_indices = MakeShapedIndices(
                centroid,
                Vector3<unsigned int>(inclusion.rows, inclusion.cols,
                                      inclusion.depth),
                layer);
            indices.insert(indices.end(), cube_indices.begin(),
                           cube_indices.end());
        }
    }

    return implicit_surface;
}

auto ImplicitSurfaceGenerator::CheckShapeUniformity() -> bool {
    const auto layers = implicit_surface.Dimension(0);

    for (int i = 0; i < layers; ++i) {
        const MatrixXr layer = implicit_surface.At(i);
        if (!CheckLayerPadding(layer)) {
            return false;
        }
    }

    return true;
}

auto ImplicitSurfaceGenerator::LayerContainsSecondaryMaterial(
    const MatrixXr& layer) -> bool {
    const auto rows = layer.rows();
    const auto cols = layer.cols();

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            if (layer(row, col) == BinaryMaterial::kSecondaryMaterial) {
                return true;
            }
        }
    }

    return false;
}

/*
@brief Ensures padding is even on all sides, returning false if not. Only
applies to isotropic material meshes with cube-like secondary material
inclusions.

This function sucks, don't have time to fix it :shrug:
*/
auto ImplicitSurfaceGenerator::CheckLayerPadding(const MatrixXr& layer)
    -> bool {
    if (!LayerContainsSecondaryMaterial(layer)) {
        return true;
    }

    const auto rows = layer.rows();
    const auto cols = layer.cols();

    int padding_amount = -1;
    int secondary_material_width = 0;

    for (int row = 0; row < rows; ++row) {
        if (padding_amount > -1) {
            break;
        }
        for (int col = 0; col < cols; ++col) {
            // First instance of the secondary material, check material width
            if (padding_amount == -1 &&
                layer(row, col) == BinaryMaterial::kSecondaryMaterial) {
                // Padding amount is left-distance on iteration.
                padding_amount = col;
                while (layer(row, col) == BinaryMaterial::kSecondaryMaterial) {
                    utils::runtime::DeformoAssert(
                        col != cols, "OVERRUN IN COLUMN RANGE: ", row, col);
                    ++secondary_material_width;
                    ++col;
                }
            }
        }
    }

    assert(padding_amount > 0 && "FAILED TO INFER PADDING");

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < padding_amount; ++col) {
            if (layer(row, col) == BinaryMaterial::kSecondaryMaterial) {
                return false;
            }
        }

        for (int col = cols - 1; col > cols - padding_amount; --col) {
            if (layer(row, col) == BinaryMaterial::kSecondaryMaterial) {
                return false;
            }
        }
    }

    return true;
}

auto ImplicitSurfaceGenerator::MakeShapedIndices(
    const Vector2<unsigned int>& centroid, const Vector3<unsigned int>& shape,
    const unsigned int layer_number) -> std::vector<Vector3<unsigned int>> {
    std::vector<Vector3<unsigned int>> indices;
    const unsigned int width = shape.x();
    const unsigned int height = shape.y();
    const unsigned int depth = shape.z();

    for (auto layer = layer_number; layer < layer_number + depth; ++layer) {
        const unsigned int current_x = centroid.x();
        const unsigned int current_y = centroid.y();
        for (auto x = 0u; x < width; ++x) {
            for (auto y = 0u; y < height; ++y) {
                indices.emplace_back(layer, current_x + x, current_y + y);
                indices.emplace_back(layer, current_x - x, current_y - y);
                indices.emplace_back(layer, current_x - x, current_y + y);
                indices.emplace_back(layer, current_x + x, current_y - y);
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

auto ImplicitSurfaceGenerator::SetFromIndices(
    const std::vector<Vector3<unsigned int>>& indices) -> void {
    for (const auto& index : indices) {
        implicit_surface(0, index.x(), index.y(), index.z());
    }
}
