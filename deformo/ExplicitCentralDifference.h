#pragma once

#include "Numerics.h"

class ExplicitCentralDifferenceMethod {
  public:
    /**
     * \brief Integration constant one, a0 = 1 / dt^2
     */
    Real a0;

    /**
     * \brief Integration constant two, a1 = 1 / 2 * dt
     */
    Real a1;

    /**
     * \brief Integration constant three, a2 = 2 * a0
     */
    Real a2;

    /**
     * \brief Integration constant four, a3 = 1 / a2
     */
    Real a3;

    /**
     * \brief Change in time per interval
     */
    const Real dt;

    /**
     * \brief The system displacement from the previous time step
     */
    VectorXr previous_position;

    ExplicitCentralDifferenceMethod(Real dt, Real point_mass,
                                    MatrixXr stiffness,
                                    const VectorXr& initial_displacements,
                                    const VectorXr& initial_forces);

    ExplicitCentralDifferenceMethod::ExplicitCentralDifferenceMethod(
        const Real dt, const SparseMatrixXr mass_matrix, MatrixXr stiffness,
        const VectorXr& initial_displacements, const VectorXr& initial_forces);

    // Setters
    void SetDamping(Real mu = 0.5f, Real lambda = 0.5f);
    void SetMassMatrix(Real point_mass);
    void SetMassMatrix(const SparseMatrixXr& m);
    void SetIntegrationConstants(Real dt) noexcept;

    // Getters
    [[nodiscard]] Real NodalMass() { return mass_matrix_.coeff(0, 0); }

    /**
    \brief Calculates the explicit Central Difference Method integration
    equation given the local position and velocity vectors.

    \n It solves according to the following algorithm:
    \n foreach time step {
        \n 1. Calculate effective loads
        \n effective_load <- forces - (stiffness - a2 * mass_matrix) *
    current_displacement - (a0 * mass_matrix) * previous_position

        \n 2. Calculate Displacements at dt
        effective_mass_matrix * effective_load = next_displacement
    \n }

    \param displacements The new displacement value
    \param forces The generalized stacked force vector
    solving for the displacement
    **/
    void Solve(VectorXr& displacements, const VectorXr& forces);

    VectorXr Velocity() const { return velocity_; }
    VectorXr Acceleration() const { return acceleration_; }

  private:
    const MatrixXr stiffness_;

    SparseMatrixXr mass_matrix_;

    MatrixXr damping_;
    SparseMatrixXr effective_mass_matrix_;

    VectorXr velocity_;
    VectorXr acceleration_;

    // Highly-specific for effective load calc
    MatrixXr el_stiffness_mass_diff_;
    MatrixXr el_mass_matrix_damping_diff_;

    void SetEffectiveMassMatrix();
    void SetLastPosition(const VectorXr& positions);
    void SetMovementVectors(const VectorXr& positions, const VectorXr& forces,
                            const MatrixXr& mass_matrix);

    VectorXr ComputeEffectiveLoad(const VectorXr& displacements,
                                  const VectorXr& forces) const;
};
