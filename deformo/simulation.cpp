#include "Simulation.h"

#include <algorithm>

Simulation::Simulation(Eigen::Ref<const Eigen::VectorXd> particles_,
                       const std::vector<Eigen::Vector3d>& fixed_indices)
    : particles(particles_) {
  // Every spring must have two nodes.
  assert(particles.size() % 2 == 0);

  AssembleFixedPointSelectionMatrices(fixed_indices);
}

void Simulation::Update() { current_time += timestep_size; }

void Simulation::Integrate() {}

void Simulation::AssembleForces() {}

void Simulation::AssembleStiffness() {}

void Simulation::AssembleMassMatrix() {
  M.resize(particles.rows(), particles.rows());
  M.setIdentity();
  M *= point_mass;
}

void Simulation::AssembleFixedPointSelectionMatrices(
    const std::vector<Eigen::Vector3d>& fixed_indices) {
  std::vector<Eigen::Triplet<double>> fixed_matrices;

  for (unsigned int i = 0; i < particles.size(); i += 3) {
    const auto vertices =
        Eigen::Vector3d(particles(i), particles(i + 1), particles(i + 2));

    // Skip over if not present.
    if (std::find(fixed_indices.begin(), fixed_indices.end(), vertices) ==
        fixed_indices.end()) {
      continue;
    }

    // Now set the triplet identity in the positions which match the vector.
    fixed_matrices.push_back(Eigen::Triplet<double>(i, i, 1.));
    fixed_matrices.push_back(Eigen::Triplet<double>(i + 1, i + 1, 1.));
    fixed_matrices.push_back(Eigen::Triplet<double>(i + 2, i + 2, 1.));
  }

  P.setFromTriplets(fixed_matrices.begin(), fixed_matrices.end());
}
