#include "Material.h"

Material MaterialFromLameCoefficients(unsigned int number,
                                      const std::string& name, Real G,
                                      Real lambda) {
    Material m;
    m.number = number;
    m.name = name;
    m.G = G;
    m.lambda = lambda;
    m.E = (G * (3 * lambda + 2 * G)) / (lambda + G);
    m.v = lambda / (2 * (lambda + G));

    return m;
}

Material MaterialFromEandv(unsigned int number, const std::string& name, Real E,
                           Real v) {
    Material m;
    m.number = number;
    m.name = name;
    m.E = E;
    m.v = v;
    m.G = E / (2 * (1 + v));
    m.lambda = ((E * v) / (1 + v) * (1 - 2 * v));

    return m;
}
