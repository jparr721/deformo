#pragma once

#include <Eigen/Dense>
#include "EigenTypes.h"

namespace integrators {
/**
@brief Calculates the explicit forward eulerian integration equation given the
local position and velocity vectors.

@param positions The generalized stacked position vector (x0, y0, z0, x1, y1, z1)
@param velocity The generalized stacked velocity vector (x0, y0, z0, x1, y1, z1)
@param forces The generalized stacked force vector (x0, y0, z0, x1, y1, z1)
@param masses The generalized mass matrix
**/
void ExplicitForwardEuler(Eigen::Vector6d& positions, Eigen::Vector6d& velocity,
                          Eigen::Ref<const Eigen::Vector6d> forces,
                          Eigen::Ref<const Eigen::Matrix66d> masses, double dt);
}  // namespace integrators
