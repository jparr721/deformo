#include "Rayleigh.h"

void ComputeRayleighDamping(MatrixXr& out, const MatrixXr& stiffness,
                            const MatrixXr& mass, Real mu, Real lambda,
                            Real modifier) {
    out = modifier * (mu * mass + lambda * stiffness);
}
