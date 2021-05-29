#pragma once

#include <Eigen/Dense>

#include "EigenTypes.h"

namespace integrators {
/**
@brief Calculates the explicit Central Difference Method integration equation
given the local position and velocity vectors.

@param positions The generalized stacked position vector
@param forces The generalized stacked force vector
@param masses The generalized triangular mass matrix
**/
void ExplicitCentralDifference(Eigen::VectorXf& displacement,
                               const Eigen::VectorXf& forces,
                               const Eigen::FullPivLU<Eigen::MatrixXf>& M_hat);
}  // namespace integrators
