#include "Simulation.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <iostream>

#include "Integrators.h"

Simulation::Simulation(double mass_, double E_, double NU_,
                       const Eigen::VectorXd& displacements_)
    : E(E_), NU(NU_), mass(mass_), displacements(displacements_) {
  InitializeIntegrationConstants();

  // Assume no initial forces so U^-dt is always the same.
  last_displacement = displacements;

  AssembleForces();
  AssembleElementStiffness();
  AssembleGlobalStiffness();
  AssembleMassMatrix();
}

void Simulation::Update() {
  current_time += timestep_size;

  // R_hat = F_ext - (K - a2M)U - (M_hat)last_displacement
  R_hat = F_ext - (K - a2 * M) * displacements - M_hat * last_displacement;

  // Set U, U', and U''
  Integrate();
}

void Simulation::Integrate() {
  const auto current_displacement = displacements;

  integrators::ExplicitCentralDifference(displacements, F_ext, M_hat_inverse);

  acceleration = a1 * (last_displacement + displacements);
  velocity =
      a0 * (last_displacement - 2 * current_displacement + displacements);

  last_displacement = current_displacement;
}

void Simulation::AssembleForces() {
  F_ext = Eigen::VectorXd::Zero(displacements.rows());
}

/**
@brief Assemble 6x6 element stiffness matrix. Given by [k] = tA[B]^T[D][B]
**/
void Simulation::AssembleElementStiffness() {
  // Clear all element stiffness values
  k.clear();

  // Node xy positions in the triangle
  double xi, yi, xj, yj, xm, ym;

  // Area
  double A;

  // Rows/Column values of the B matrix
  double betai, betaj, betam, gammai, gammaj, gammam;

  // Plane interaction matrix
  Eigen::Matrix36d B;

  // Plane stress/strain matrix
  Eigen::Matrix3d D;

  // Iterate by groups of 3.
  for (std::size_t i = 0; i < displacements.size(); i += 6) {
    xi = displacements[i];
    yi = displacements[i + 1];
    xj = displacements[i + 2];
    yj = displacements[i + 3];
    xm = displacements[i + 4];
    ym = displacements[i + 5];

    A = (xi * (yj - ym) + xj * (ym - yi) + xm * (yi - yj)) / 2;

    betai = yj = ym;
    betaj = ym - yi;
    betam = yi - yj;

    gammai = xm - xj;
    gammaj = xi - xm;
    gammam = xj - xi;

    // Construct beta as a 6x6 matrix.
    B << betai, 0, betaj, 0, betam, 0;
    B << 0, gammai, 0, gammaj, 0, gammam;
    B << gammai, betai, gammaj, betaj, gammam, betam;
    B /= (2 * A);

    // Plane stress calculation
    D << 1, NU, 0;
    D << NU, 1, 0;
    D << 0, 0, (1 - NU) / 2;

    // Plane stress calculation
    D *= E / (1 - std::pow(NU, 2));

    const uint64_t node_number = i / 3;
    const Eigen::Matrix66d kk = t * A * B.transpose() * D * B;
    k.emplace_back(std::make_pair(
        kk,
        std::vector<uint64_t>{node_number, node_number + 1, node_number + 2}));
  }
}

void Simulation::AssembleGlobalStiffness() {
  const double rowsize = displacements.rows() * 2;

  // Compute row size with 3 degrees of freedom.
  K.resize(rowsize, rowsize);

  // Iterate through all points and element stiffness matrices, generating the
  // output inside of the K matrix.
  for (const auto& [kk, pts] : k) {
    const auto i = pts[0];
    const auto j = pts[1];
    const auto m = pts[2];

    // oof.
    K(2 * i - 1, 2 * i - 1) = K(2 * i - 1, 2 * i - 1) + kk(1, 1);
    K(2 * i - 1, 2 * i) = K(2 * i - 1, 2 * i) + kk(1, 2);
    K(2 * i - 1, 2 * j - 1) = K(2 * i - 1, 2 * j - 1) + kk(1, 3);
    K(2 * i - 1, 2 * j) = K(2 * i - 1, 2 * j) + kk(1, 4);
    K(2 * i - 1, 2 * m - 1) = K(2 * i - 1, 2 * m - 1) + kk(1, 5);
    K(2 * i - 1, 2 * m) = K(2 * i - 1, 2 * m) + kk(1, 6);
    K(2 * i, 2 * i - 1) = K(2 * i, 2 * i - 1) + kk(2, 1);
    K(2 * i, 2 * i) = K(2 * i, 2 * i) + kk(2, 2);
    K(2 * i, 2 * j - 1) = K(2 * i, 2 * j - 1) + kk(2, 3);
    K(2 * i, 2 * j) = K(2 * i, 2 * j) + kk(2, 4);
    K(2 * i, 2 * m - 1) = K(2 * i, 2 * m - 1) + kk(2, 5);
    K(2 * i, 2 * m) = K(2 * i, 2 * m) + kk(2, 6);
    K(2 * j - 1, 2 * i - 1) = K(2 * j - 1, 2 * i - 1) + kk(3, 1);
    K(2 * j - 1, 2 * i) = K(2 * j - 1, 2 * i) + kk(3, 2);
    K(2 * j - 1, 2 * j - 1) = K(2 * j - 1, 2 * j - 1) + kk(3, 3);
    K(2 * j - 1, 2 * j) = K(2 * j - 1, 2 * j) + kk(3, 4);
    K(2 * j - 1, 2 * m - 1) = K(2 * j - 1, 2 * m - 1) + kk(3, 5);
    K(2 * j - 1, 2 * m) = K(2 * j - 1, 2 * m) + kk(3, 6);
    K(2 * j, 2 * i - 1) = K(2 * j, 2 * i - 1) + kk(4, 1);
    K(2 * j, 2 * i) = K(2 * j, 2 * i) + kk(4, 2);
    K(2 * j, 2 * j - 1) = K(2 * j, 2 * j - 1) + kk(4, 3);
    K(2 * j, 2 * j) = K(2 * j, 2 * j) + kk(4, 4);
    K(2 * j, 2 * m - 1) = K(2 * j, 2 * m - 1) + kk(4, 5);
    K(2 * j, 2 * m) = K(2 * j, 2 * m) + kk(4, 6);
    K(2 * m - 1, 2 * i - 1) = K(2 * m - 1, 2 * i - 1) + kk(5, 1);
    K(2 * m - 1, 2 * i) = K(2 * m - 1, 2 * i) + kk(5, 2);
    K(2 * m - 1, 2 * j - 1) = K(2 * m - 1, 2 * j - 1) + kk(5, 3);
    K(2 * m - 1, 2 * j) = K(2 * m - 1, 2 * j) + kk(5, 4);
    K(2 * m - 1, 2 * m - 1) = K(2 * m - 1, 2 * m - 1) + kk(5, 5);
    K(2 * m - 1, 2 * m) = K(2 * m - 1, 2 * m) + kk(5, 6);
    K(2 * m, 2 * i - 1) = K(2 * m, 2 * i - 1) + kk(6, 1);
    K(2 * m, 2 * i) = K(2 * m, 2 * i) + kk(6, 2);
    K(2 * m, 2 * j - 1) = K(2 * m, 2 * j - 1) + kk(6, 3);
    K(2 * m, 2 * j) = K(2 * m, 2 * j) + kk(6, 4);
    K(2 * m, 2 * m - 1) = K(2 * m, 2 * m - 1) + kk(6, 5);
    K(2 * m, 2 * m) = K(2 * m, 2 * m) + kk(6, 6);
  }
}

void Simulation::InitializeVelocity() {
  velocity = Eigen::VectorXd::Zero(displacements.rows());
}

void Simulation::InitializeAcceleration() {
  acceleration = Eigen::VectorXd::Zero(displacements.rows());
}

void Simulation::AssembleMassMatrix() {
  std::vector<Eigen::Triplet<double>> mass_entries;
  mass_entries.reserve(displacements.rows());

  M.resize(displacements.rows(), displacements.rows());

  for (int i = 0; i < displacements.rows(); ++i) {
    mass_entries.push_back(Eigen::Triplet<double>(i, i, mass));
  }

  M.setFromTriplets(mass_entries.begin(), mass_entries.end());

  // Set effective mass matrix
  M_hat = (a0 * M).eval();

  // Triangularize M_hat
  Eigen::SimplicialLDLT<Eigen::SparseMatrixXd> solver;
  solver.compute(M_hat);

  // If the solver fails for some reason, regularize the result and re-run
  if (solver.info() != Eigen::Success) {
    std::cout << "Failed to solve M_hat" << std::endl;
  }

  M_hat_inverse = solver.solve(M_hat);
}

void Simulation::InitializeIntegrationConstants() {
  a0 = 1 / std::pow(timestep_size, 2);
  a1 = 1 / (timestep_size * 2);
  a2 = 2 * a0;
  a3 = 1 / a2;
}
