#include "Integrators.h"
#include "Utils.h"

void ExplicitCentralDifferenceMethod::Solve(Eigen::VectorXf& positions,
                                            const Eigen::VectorXf& forces) {
    const Eigen::VectorXf effective_load =
        forces - (stiffness_ - a2 * mass_matrix_) * positions -
        (a0 * mass_matrix_)*previous_position;

    const Eigen::VectorXf next_displacement =
        effective_mass_matrix_.inverse() * effective_load;

    acceleration_ =
        a0 * (previous_position - 2 * positions + next_displacement);
    velocity_ = a1 * ((-1 * previous_position) + next_displacement);
    previous_position = positions;
    positions = next_displacement;
}

void ExplicitCentralDifferenceMethod::SetLastPosition(
    const Eigen::VectorXf& positions) {
    previous_position = positions - dt * velocity_ + a3 * acceleration_;
}
