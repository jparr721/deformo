#include "BoundaryCondition.h"

BoundaryConditions AssignBoundaryConditionToFixedNodes(
    const std::vector<unsigned int>& face_indices,
    const Eigen::Vector3f& force) {
    std::vector<BoundaryCondition> boundary_conditions;
    for (const auto& face : face_indices) {
        boundary_conditions.emplace_back(BoundaryCondition{face, force});
    }

    return boundary_conditions;
}