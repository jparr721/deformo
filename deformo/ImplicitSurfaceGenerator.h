#pragma once

#include "Numerics.h"
#include <vector>

struct BinaryInclusion {
    unsigned int n_inclusions;
    unsigned int area;
    unsigned int depth;
    unsigned int rows;
    unsigned int cols;
};

class ImplicitSurfaceGenerator {
  public:
    enum BinaryMaterial {
        kSecondaryMaterial = 0x00,
        kPrimaryMaterial = 0x01,
    };

    int minimum_surface_padding = 1;

    Tensor3r implicit_surface;

    ImplicitSurfaceGenerator(const unsigned int height,
                             const unsigned int width,
                             const unsigned int depth);
    auto Generate(const BinaryInclusion inclusion) -> Tensor3r;

  private:
    // Different Isotropic Material Checks
    auto CheckShapeUniformity() -> bool;

    // Helpers
    auto LayerContainsSecondaryMaterial(const MatrixXr& layer) -> bool;
    auto CheckLayerPadding(const MatrixXr& layer) -> bool;

    // Shape Generators
    auto MakeShapedIndices(const Vector2<unsigned>& centroid,
                           const Vector3<unsigned>& shape,
                           unsigned int layer_number)
        -> std::vector<Vector3<unsigned int>>;
    auto SetFromIndices(const std::vector<Vector3<unsigned int>>& indices) -> void;
};
