#pragma once

#include "Numerics.h"
#include <string>

struct Material {
    // The number this material represents in the scalar field
    unsigned int number = 1;

    // Name of the material, for bookkeeping
    std::string name = "";

    // Young's Modulus
    Real E = -1;

    // Poisson's Ratio
    Real v = -1;

    // Shear Modulus
    Real G = -1;

    // Lame's Lambda
    Real lambda = -1;

    auto IsInit() const noexcept -> bool { return !name.empty(); }
};

Material MaterialFromLameCoefficients(unsigned int number,
                                      const std::string& name, Real G,
                                      Real lambda);

Material MaterialFromEandv(unsigned int number, const std::string& name, Real E,
                           Real v);

