#pragma once

#include "EigenTypes.h"
#include <utility>

class ExplicitCentralDifferenceMethod {
  public:
    /**
     * \brief Integration constant one, a0 = 1 / dt^2
     */
    float a0;

    /**
     * \brief Integration constant two, a1 = 1 / 2 * dt
     */
    float a1;

    /**
     * \brief Integration constant three, a2 = 2 * a0
     */
    float a2;

    /**
     * \brief Integration constant four, a3 = 1 / a2
     */
    float a3;

    /**
     * \brief Change in time per interval
     */
    const float dt;

    /**
     * \brief The system displacement from the previous time step
     */
    Eigen::VectorXf previous_position;

    // TODO(@jpar721) - Add damping to effective matrix calc.
    ExplicitCentralDifferenceMethod(const float dt,
                                    const Eigen::VectorXf& displacements,
                                    Eigen::MatrixXf stiffness,
                                    const Eigen::SparseMatrixXf& mass_matrix,
                                    const Eigen::VectorXf& initial_forces,
                                    const Eigen::MatrixXf& damping);

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
    void Solve(Eigen::VectorXf& displacements, const Eigen::VectorXf& forces);

    const Eigen::VectorXf Velocity() const { return velocity_; }
    const Eigen::VectorXf Acceleration() const { return acceleration_; }

  private:
    const Eigen::MatrixXf damping_;
    const Eigen::MatrixXf stiffness_;
    const Eigen::SparseMatrixXf mass_matrix_;
    Eigen::SparseMatrixXf effective_mass_matrix_;

    Eigen::VectorXf velocity_;
    Eigen::VectorXf acceleration_;

    // Highly-specific for effective load calc
    Eigen::MatrixXf el_stiffness_mass_diff_;
    Eigen::MatrixXf el_mass_matrix_damping_diff_;

    void SetLastPosition(const Eigen::VectorXf& positions);
    void SetEffectiveMassMatrix();
    void SetEffectiveLoadConstants();
    void SetIntegrationConstants();
    void SetMovementVectors(const Eigen::VectorXf& positions,
                            const Eigen::VectorXf& forces,
                            const Eigen::MatrixXf& mass_matrix);

    Eigen::VectorXf ComputeEffectiveLoad(const Eigen::VectorXf& displacements,
                                         const Eigen::VectorXf& forces);
};
