#include "ExplicitCentralDifference.h"
#include "Rayleigh.h"
#include "Utils.h"

ExplicitCentralDifferenceMethod::ExplicitCentralDifferenceMethod(
    const Real dt, const SparseMatrixXr mass_matrix, MatrixXr stiffness,
    const VectorXr& initial_displacements, const VectorXr& initial_forces)
    : dt(dt), stiffness_(std::move(stiffness)) {
    SetMassMatrix(mass_matrix);
    SetDamping(0.5f, 0.5f);
    SetIntegrationConstants(dt);
    SetEffectiveMassMatrix();
    SetMovementVectors(initial_displacements, initial_forces, mass_matrix_);
    SetLastPosition(initial_displacements);
}

ExplicitCentralDifferenceMethod::ExplicitCentralDifferenceMethod(
    const Real dt, const Real point_mass, MatrixXr stiffness,
    const VectorXr& initial_displacements, const VectorXr& initial_forces)
    : dt(dt), stiffness_(std::move(stiffness)) {
    SetMassMatrix(point_mass);
    SetDamping(0.5f, 0.5f);
    SetIntegrationConstants(dt);
    SetEffectiveMassMatrix();
    SetMovementVectors(initial_displacements, initial_forces, mass_matrix_);
    SetLastPosition(initial_displacements);
}

void ExplicitCentralDifferenceMethod::SetDamping(const Real mu,
                                                 const Real lambda) {
    ComputeRayleighDamping(damping_, stiffness_, mass_matrix_, mu, lambda, 0);
}

void ExplicitCentralDifferenceMethod::Solve(VectorXr& displacements,
                                            const VectorXr& forces) {
    const VectorXr effective_load = ComputeEffectiveLoad(displacements, forces);

    const VectorXr next_displacement = effective_mass_matrix_ * effective_load;

    acceleration_ =
        a0 * (previous_position - 2 * displacements + next_displacement);
    velocity_ = a1 * ((-1 * previous_position) + next_displacement);
    previous_position = displacements;
    displacements = next_displacement;
}

void ExplicitCentralDifferenceMethod::SetLastPosition(
    const VectorXr& positions) {
    previous_position = positions - dt * velocity_ + a3 * acceleration_;
}

void ExplicitCentralDifferenceMethod::SetMassMatrix(Real point_mass) {
    mass_matrix_.resize(stiffness_.rows(), stiffness_.cols());
    mass_matrix_.setIdentity();
    mass_matrix_ *= point_mass;
}

void ExplicitCentralDifferenceMethod::SetMassMatrix(const SparseMatrixXr& m) {
    mass_matrix_ = m;
}

void ExplicitCentralDifferenceMethod::SetEffectiveMassMatrix() {
    effective_mass_matrix_ = a0 * mass_matrix_ + a1 * damping_;
    Eigen::SparseLU<SparseMatrixXr> solver;
    solver.compute(effective_mass_matrix_);
    SparseMatrixXr I(effective_mass_matrix_.rows(),
                     effective_mass_matrix_.cols());
    I.setIdentity();
    effective_mass_matrix_ = solver.solve(I);
}

void ExplicitCentralDifferenceMethod::SetIntegrationConstants(
    const Real dt) noexcept {
    a0 = 1.f / (std::powf(dt, 2));
    a1 = 1.f / (2.f * dt);
    a2 = 2.f * a0;
    a3 = 1.f / a2;
}

void ExplicitCentralDifferenceMethod::SetMovementVectors(
    const VectorXr& positions, const VectorXr& forces,
    const MatrixXr& mass_matrix) {
    velocity_.resize(positions.rows());
    velocity_.setZero();
    acceleration_.resize(positions.rows());
    acceleration_ = mass_matrix.inverse() * forces;
}

VectorXr ExplicitCentralDifferenceMethod::ComputeEffectiveLoad(
    const VectorXr& displacements, const VectorXr& forces) const {
    return forces - (stiffness_ - a2 * mass_matrix_) * displacements -
           (a0 * mass_matrix_ - a1 * damping_) * previous_position;
}
