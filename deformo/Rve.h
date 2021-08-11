#pragma once

#include "ImplicitSurfaceGenerator.h"
#include "Numerics.h"
#include "Material.h"
#include <string>
#include <utility>

class Rve {
  public:
    // Homogenous 1-material structure
    const bool homogenous;

    // RVE Dimenions
    const unsigned int width;
    const unsigned int height;
    const unsigned int depth;

    // RVE Inclusion Dimensions
    unsigned int inclusion_width;
    unsigned int inclusion_height;
    unsigned int inclusion_depth;

    // The applied strain
    Real strain;

    // Tetgen tetrahedral element volume
    Real mesh_density;

    // The fraction of the RVE which contains material 2
    Real volume_fraction;

    Material material_1;
    Material material_2;

    Rve(const bool homogenous, const unsigned int width,
        const unsigned int height, const unsigned int depth, const Real strain,
        const Real mesh_density, const Real volume_fraction)
        : homogenous(homogenous), width(width), height(height), depth(depth),
          strain(strain), mesh_density(mesh_density),
          volume_fraction(volume_fraction) {}

    auto ToImplicitSurface() -> Tensor3r;
};
