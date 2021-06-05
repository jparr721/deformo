#pragma once

#include "EigenTypes.h"
#include <Eigen/Sparse>

class ExplicitCentralDifferenceMethod {
  public:
    /**
     * \brief Integration constant one, a0 = 1 / dt^2
     */
    const float a0;

    /**
     * \brief Integration constant two, a1 = 1 / 2 * dt
     */
    const float a1;

    /**
     * \brief Integration constant three, a2 = 2 * a0
     */
    const float a2;

    /**
     * \brief The system displacement from the previous timestep
     */
    Eigen::VectorXf previous_position;

    ExplicitCentralDifferenceMethod(const float a0, const float a1,
                                    const float a2,
                                    const Eigen::VectorXf& initial_displacement,
                                    const Eigen::MatrixXf& stiffness,
                                    const Eigen::SparseMatrixXf& mass_matrix)
        : a0(a0), a1(a1), a2(a2), previous_position(initial_displacement),
          stiffness_(stiffness), mass_matrix_(mass_matrix),
          effective_mass_matrix_(
              Eigen::SimplicialLDLT(Eigen::SparseMatrixXf(a0 * mass_matrix))) {}

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
    \param acceleration The acceleration of the positions at timestep t
    \param velocity The velocity of the positions at timestep t
    \param forces The generalized stacked force vector
    solving for the displacement
    **/
    void Solve(Eigen::VectorXf& positions, Eigen::VectorXf& acceleration,
               Eigen::VectorXf& velocity, const Eigen::VectorXf& forces);

  private:
    const Eigen::MatrixXf stiffness_;
    const Eigen::SparseMatrixXf mass_matrix_;
    const Eigen::SimplicialLDLT<Eigen::SparseMatrixXf> effective_mass_matrix_;
};
