#pragma once

#include "EigenTypes.h"
#include <Eigen/Sparse>
#include <igl/barycenter.h>

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
                                    const Eigen::VectorXf& positions,
                                    const Eigen::MatrixXf& stiffness,
                                    const Eigen::SparseMatrixXf& mass_matrix)
        : dt(dt), stiffness_(stiffness), mass_matrix_(mass_matrix) {
        a0 = 1.f / (std::powf(dt, 2));
        a1 = 1.f / (2.f * dt);
        a2 = 2.f * a0;
        a3 = 1.f / a2;

        velocity_.resize(positions.rows());
        velocity_.setZero();
        acceleration_.resize(positions.rows());
        acceleration_.setZero();
        effective_mass_matrix_ = a0 * mass_matrix_;
        SetLastPosition(positions);
    }

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

    \param positions The new displacement value
    \param forces The generalized stacked force vector
    solving for the displacement
    **/
    void Solve(Eigen::VectorXf& positions, const Eigen::VectorXf& forces);

  private:
    const Eigen::MatrixXf stiffness_;
    const Eigen::SparseMatrixXf mass_matrix_;
    Eigen::MatrixXf effective_mass_matrix_;

    Eigen::VectorXf velocity_;
    Eigen::VectorXf acceleration_;

    void SetLastPosition(const Eigen::VectorXf& positions);
};
