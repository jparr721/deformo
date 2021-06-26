#pragma once

#include <Eigen/Dense>
#include <vector>

/*
Static boundary condition for a given node
*/
struct BoundaryCondition {
    unsigned int node;
    Eigen::Vector3f force;
};

using BoundaryConditions = std::vector<BoundaryCondition>;

BoundaryConditions AssignBoundaryConditionToFixedNodes(
    const std::vector<unsigned int>& face_indices,
    const Eigen::Vector3f& force);
