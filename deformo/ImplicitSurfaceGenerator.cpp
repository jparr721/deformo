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

ImplicitSurfaceGenerator::ImplicitSurfaceGenerator(unsigned int height,
                                                   unsigned int width,
                                                   unsigned int depth) {
    implicit_surface.Resize(height, width, depth);
    implicit_surface.Instance().setConstant(BinaryMaterial::kPrimaryMaterial);
}

auto ImplicitSurfaceGenerator::Generate(const BinaryInclusion inclusion)
    -> Tensor3r {
    // Y axis origin with padding
    VectorXr y_axis_origins =
        linear_algebra::LinSpace(inclusion.area + minimum_surface_padding,
                                 (implicit_surface.Dimension(2) + 1 -
                                  (inclusion.area + minimum_surface_padding)),
                                 inclusion.rows);
    y_axis_origins -= VectorXr::Ones(y_axis_origins.rows());
    CheckLayerPadding(y_axis_origins, implicit_surface.Dimension(2));

    VectorXr x_axis_origins =
        linear_algebra::LinSpace(inclusion.area + minimum_surface_padding,
                                 (implicit_surface.Dimension(1) + 1 -
                                  (inclusion.area + minimum_surface_padding)),
                                 inclusion.cols);
    x_axis_origins -= VectorXr::Ones(x_axis_origins.rows());
    CheckLayerPadding(x_axis_origins, implicit_surface.Dimension(1));

    std::vector<Vector2<unsigned int>> centroids;
    for (int i = 0; i < y_axis_origins.rows(); ++i) {
        for (int j = 0; j < x_axis_origins.rows(); ++j) {
            const int y = static_cast<int>(y_axis_origins(i));
            const int x = static_cast<int>(x_axis_origins(j));
            centroids.emplace_back(x, y);
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
                Vector3<unsigned int>(inclusion.area, inclusion.area,
                                      inclusion.depth),
                layer);
            indices.insert(indices.end(), cube_indices.begin(),
                           cube_indices.end());
        }
    }

    SetFromIndices(indices);

    return implicit_surface;
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
applies to isotropic material meshes.
*/
auto ImplicitSurfaceGenerator::CheckLayerPadding(const VectorXr& layer,
                                                 const int max) -> void {
    const auto left_padding = layer(0) + 1;
    const auto right_padding = max - layer(layer.rows() - 1);

    deformo_assert.Assert(left_padding == right_padding, __FUNCTION__, __FILE__,
                          __LINE__, "Padding does not match: ", left_padding,
                          right_padding, layer);
}

auto ImplicitSurfaceGenerator::MakeShapedIndices(
    const Vector2<unsigned int>& centroid, const Vector3<unsigned int>& shape,
    const unsigned int layer_number) -> std::vector<Vector3<unsigned int>> {
    std::vector<Vector3<unsigned int>> indices;
    const unsigned int width = static_cast<unsigned int>(shape.x() / 2);
    const unsigned int height = static_cast<unsigned int>(shape.y() / 2);
    const unsigned int depth = shape.z();

    const unsigned int current_x = centroid.x();
    const unsigned int current_y = centroid.y();
    for (auto layer = layer_number; layer < layer_number + depth; ++layer) {
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
    for (const Vector3<unsigned int>& index : indices) {
        implicit_surface(BinaryMaterial::kSecondaryMaterial, index.x(),
                         index.y(), index.z());
    }
}
