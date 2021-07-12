#pragma once

#include "Numerics.h"

class ImplicitSurfaceGenerator {
  public:
    const unsigned int height;
    const unsigned int width;
    const unsigned int depth;

    Tensor3r implicit_surface;

    ImplicitSurfaceGenerator(const unsigned int height, const unsigned int width,
                             const unsigned int depth)
        : height(height), width(width), depth(depth) {}
    auto Generate() -> Tensor3r;

private:
    auto ValidatePostGenerationChecks() -> bool;

    // Different Isotropic Material Checks
    auto CheckShapeUniformity() -> bool;
    auto CheckShapeConsistency() -> bool;

    // Helpers
    auto LayerContainsVoid(const MatrixXr& layer) -> bool;
    auto CheckLayerPadding(const MatrixXr& layer) -> bool;
    auto CheckLayerToLayerPadding() -> bool;
};
