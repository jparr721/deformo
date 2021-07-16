#pragma once

#include "Numerics.h"
#include <string>
#include <utility>

struct Material {
    // Name of the material, for bookkeeping
    std::string name = "";

    // Young's Modulus
    Real E = -1;

    // Poisson's Ratio
    Real v = -1;

    auto Empty() -> bool;
};

class Rve {
  public:
    Real volume_fraction;
    Real strain;
    Real mesh_density;

    Material material_1;
    Material material_2;

    Rve(const Real strain, const Real mesh_density)
        : strain(strain), mesh_density(mesh_density) {}

    auto AssignSections(const std::string& material_name) -> void;
    auto SetVolumeFractionForMaterial(const std::string& material_name,
                                      Real volume_fraction) -> void;
};
