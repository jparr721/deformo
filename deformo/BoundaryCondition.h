#pragma once

#include <Eigen/Dense>

/*
Static boundary condition for a given node
*/
struct BoundaryCondition {
  unsigned int node;
  Eigen::Vector3f force;
};

