#include "Rayleigh.h"

void ComputeRayleighDamping(Eigen::MatrixXf& out,
                           const Eigen::MatrixXf& stiffness,
                           const Eigen::MatrixXf& mass, float mu, float lambda,
                           float modifier) {
    out = modifier * (mu * mass + lambda * stiffness);
}
