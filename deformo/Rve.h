#pragma once

#include "ImplicitSurfaceGenerator.h"
#include "Numerics.h"
#include <string>
#include <utility>

struct Material {
    // The number this material represents in the scalar field
    unsigned int number;

    // Name of the material, for bookkeeping
    std::string name;

    // Young's Modulus
    Real E = -1;

    // Poisson's Ratio
    Real v = -1;

    auto Lambda() -> Real {
        // TODO(@jparr721) - Implement
        return E;
    }
    auto Mu() -> Real {
        // TODO(@jparr721) - Implement
        return v;
    }
};

class Rve {
  public:
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

    Rve(const unsigned int width, const unsigned int height,
        const unsigned int depth, const Real strain, const Real mesh_density,
        const Real volume_fraction)
        : width(width), height(height), depth(depth), strain(strain),
          mesh_density(mesh_density), volume_fraction(volume_fraction) {}
    auto ToImplicitSurface() -> Tensor3r;

    /*
    @brief Set the dimensions for material_2. This might change the volume
    fraction if it cannot be reasonably normalized. These values will be
    automatically set if unset during the inclusion generation process.
    */
    auto SetInclusionDimenions(const Vector3<unsigned int>& dimensions) -> void;

  private:
    auto MakeVolumeAwareBinaryInclusion() -> BinaryInclusion;
};
