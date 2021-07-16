#include "ImplicitSurfaceGenerator.h"
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

Tensor3r ImplicitSurfaceGenerator::Generate() { return implicit_surface; }

bool ImplicitSurfaceGenerator::ValidatePostGenerationChecks() { return true; }

bool ImplicitSurfaceGenerator::CheckShapeUniformity() { return true; }

bool ImplicitSurfaceGenerator::CheckShapeConsistency() { return true; }

bool ImplicitSurfaceGenerator::LayerContainsVoid(const MatrixXr& layer) {
    return true;
}

bool ImplicitSurfaceGenerator::CheckLayerPadding(const MatrixXr& layer) {
    return true;
}

bool ImplicitSurfaceGenerator::CheckLayerToLayerPadding() { return true; }

std::vector<Vector3<unsigned int>> ImplicitSurfaceGenerator::MakeShapedIndices(
    const Vector2<unsigned int>& centroid, const Vector3<unsigned int>& shape,
    const unsigned int layer_number) {
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
    const std::vector<Vector3<unsigned int>>& indices) {
    for (const auto& index : indices) {
        implicit_surface(0, index.x(), index.y(), index.z());
    }
}
