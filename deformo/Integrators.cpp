#include "Integrators.h"
#include "Utils.h"

void ExplicitCentralDifferenceMethod::Solve(Eigen::VectorXf& positions,
                                            Eigen::VectorXf& acceleration,
                                            Eigen::VectorXf& velocity,
                                            const Eigen::VectorXf& forces) {
    const Eigen::VectorXf current_positions = positions;
    const Eigen::VectorXf effective_load =
        forces - (((stiffness_ - a2 * mass_matrix_) * positions) -
                  ((a0 * mass_matrix_) * previous_position));

    utils::GTestDebugPrint(effective_load);
    std::cerr << "==" << std::endl;
    const Eigen::VectorXf next_displacement =
        effective_mass_matrix_.solve(effective_load);

    utils::GTestDebugPrint(next_displacement);
    std::cerr << "==" << std::endl;
    positions += next_displacement;
    utils::GTestDebugPrint(positions);
    acceleration =
        a0 * (previous_position - (2 * current_positions) + positions);
    velocity = a1 * ((-1 * previous_position) + positions);
    previous_position = current_positions;
}
