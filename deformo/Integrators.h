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
void ExplicitCentralDifference(
    Eigen::VectorXd& displacement, Eigen::Ref<const Eigen::VectorXd> forces,
    Eigen::Ref<const Eigen::MatrixXd> mass_triangular);
}  // namespace integrators
