#include "Integrators.h"
#include "Utils.h"

void ExplicitCentralDifferenceMethod::Solve(Eigen::VectorXf& displacements,
                                            const Eigen::VectorXf& forces) {
    const Eigen::VectorXf effective_load =
        ComputeEffectiveLoad(displacements, forces);

    const Eigen::VectorXf next_displacement =
        effective_mass_matrix_ * effective_load;

    acceleration_ =
        a0 * (previous_position - 2 * displacements + next_displacement);
    velocity_ = a1 * ((-1 * previous_position) + next_displacement);
    previous_position = displacements;
    displacements = next_displacement;
}

void ExplicitCentralDifferenceMethod::SetLastPosition(
    const Eigen::VectorXf& positions) {
    previous_position = positions - dt * velocity_ + a3 * acceleration_;
}

void ExplicitCentralDifferenceMethod::SetEffectiveMassMatrix() {
    effective_mass_matrix_ = a0 * mass_matrix_;
    Eigen::SparseLU<Eigen::SparseMatrixXf> solver;
    solver.compute(effective_mass_matrix_);
    Eigen::SparseMatrixXf I(effective_mass_matrix_.rows(),
                            effective_mass_matrix_.cols());
    I.setIdentity();
    effective_mass_matrix_ = solver.solve(I);
}

void ExplicitCentralDifferenceMethod::SetEffectiveLoadConstants() {
    el_stiffness_mass_diff_ = stiffness_ - a2 * mass_matrix_;
    el_mass_matrix_damping_diff_ = a0 * mass_matrix_;
}

void ExplicitCentralDifferenceMethod::SetIntegrationConstants() {
    a0 = 1.f / (std::powf(dt, 2));
    a1 = 1.f / (2.f * dt);
    a2 = 2.f * a0;
    a3 = 1.f / a2;
}

void ExplicitCentralDifferenceMethod::SetMovementVectors(
    const Eigen::VectorXf& positions) {
    velocity_.resize(positions.rows());
    velocity_.setZero();
    acceleration_.resize(positions.rows());
    acceleration_.setZero();
}

Eigen::VectorXf ExplicitCentralDifferenceMethod::ComputeEffectiveLoad(
    const Eigen::VectorXf& displacements, const Eigen::VectorXf& forces) {
    return forces - el_stiffness_mass_diff_ * displacements -
           el_mass_matrix_damping_diff_ * previous_position;
}
