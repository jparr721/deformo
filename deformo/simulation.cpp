#include "Simulation.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <iostream>

#include "Integrators.h"

Simulation::Simulation(double mass_, double E_, double NU_,
                       const std::shared_ptr<Mesh>& mesh_)
    : E(E_), NU(NU_), mass(mass_), mesh(std::move(mesh_)) {
  // Assert minimum size is at least 3 so we can do triangle calculations
  assert(mesh->rows() >= 6);
  InitializeIntegrationConstants();

  const double K_rowsize = mesh->rows() * 2;

  K = Eigen::MatrixXd::Zero(K_rowsize, K_rowsize);

  // Assume no initial forces so U^-dt is always the same.
  last_displacement = mesh->vertices;

  AssembleForces();
  AssembleElementStiffness();
  AssembleGlobalStiffness();
  AssembleMassMatrix();
}

void Simulation::Update() {
  current_time += timestep_size;

  // R_hat = F_ext - (K - a2M)U - (M_hat)last_displacement
  R_hat = F_ext - (K - a2 * M) * mesh->vertices - M_hat * last_displacement;

  // Set U, U', and U''
  Integrate();
}

void Simulation::Integrate() {
  const auto current_displacement = mesh->vertices;

  integrators::ExplicitCentralDifference(mesh->vertices, F_ext, M_hat_inverse);

  acceleration = a1 * (last_displacement + mesh->vertices);
  velocity =
      a0 * (last_displacement - 2 * current_displacement + mesh->vertices);

  last_displacement = current_displacement;
}

void Simulation::AssembleForces() {
  F_ext = Eigen::VectorXd::Zero(mesh->rows());
}

/**
@brief Assemble 6x6 element stiffness matrix. Given by [k] = tA[B]^T[D][B]
**/
void Simulation::AssembleElementStiffness() {
  // Clear all element stiffness values
  k.clear();

  // Iterate by groups of 3.
  for (std::size_t i = 0; i < mesh->rows(); i += 6) {
    const double xi = mesh->vertices(i);
    const double yi = mesh->vertices(i + 1);
    const double xj = mesh->vertices(i + 2);
    const double yj = mesh->vertices(i + 3);
    const double xm = mesh->vertices(i + 4);
    const double ym = mesh->vertices(i + 5);

    const double A = (xi * (yj - ym) + xj * (ym - yi) + xm * (yi - yj)) / 2.;

    const double gammai = xm - xj;
    const double gammaj = xi - xm;
    const double gammam = xj - xi;

    const double betai = yj - ym;
    const double betaj = ym - yi;
    const double betam = yi - yj;

    // Construct beta as a 3x6 matrix.
    Eigen::Matrix36d B;
    B.row(0) << betai, 0., betaj, 0., betam, 0.;
    B.row(1) << 0., gammai, 0., gammaj, 0., gammam;
    B.row(2) << gammai, betai, gammaj, betaj, gammam, betam;
    B /= (2. * A);

    // Plane stress/strain matrix
    Eigen::Matrix3d D;
    // Plane stress calculation
    D.row(0) << 1., NU, 0.;
    D.row(1) << NU, 1, 0.;
    D.row(2) << 0., 0., (1. - NU) / 2.;

    // Plane stress calculation
    D *= E / (1. - std::pow(NU, 2));

    const uint64_t node_number = (i / 3) + 1;
    const Eigen::Matrix66d kk = t * A * B.transpose() * D * B;

    // Nodes ordered in counter-clockwise fashion
    const auto l_node_number = mesh->index(xi, yi);
    const auto m_node_number = mesh->index(xj, yj);
    const auto r_node_number = mesh->index(xm, ym);

    const ElementStiffness stiffness_entry = {
        .stiffness_matrix = kk,
        .indices = std::vector<unsigned int>{l_node_number, m_node_number,
                                             r_node_number},
    };

    k.emplace_back(stiffness_entry);
  }
}

void Simulation::AssembleGlobalStiffness() {
  // Iterate through all points and element stiffness matrices, generating the
  // output inside of the K matrix.
  for (const auto& element_stiffness : k) {
    const auto kk = element_stiffness.stiffness_matrix;
    const auto pts = element_stiffness.indices;

    const auto i = pts[0];
    const auto j = pts[1];
    const auto m = pts[2];

    // oof.
    K(2 * i - 2, 2 * i - 2) += kk(0, 0);
    K(2 * i - 2, 2 * i - 1) += kk(0, 1);
    K(2 * i - 2, 2 * j - 1) += kk(0, 2);
    K(2 * i - 2, 2 * j - 1) += kk(0, 3);
    K(2 * i - 2, 2 * m - 1) += kk(0, 4);
    K(2 * i - 2, 2 * m - 1) += kk(0, 5);
    K(2 * i - 1, 2 * i - 2) += kk(1, 0);
    K(2 * i - 1, 2 * i - 1) += kk(1, 1);
    K(2 * i - 1, 2 * j - 2) += kk(1, 2);
    K(2 * i - 1, 2 * j - 1) += kk(1, 3);
    K(2 * i - 1, 2 * m - 1) += kk(1, 4);
    K(2 * i - 2, 2 * m - 1) += kk(1, 5);
    K(2 * j - 2, 2 * i - 2) += kk(2, 0);
    K(2 * j - 2, 2 * i - 1) += kk(2, 1);
    K(2 * j - 2, 2 * j - 2) += kk(2, 2);
    K(2 * j - 2, 2 * j - 1) += kk(2, 3);
    K(2 * j - 2, 2 * m - 2) += kk(2, 4);
    K(2 * j - 2, 2 * m - 1) += kk(2, 5);
    K(2 * j - 1, 2 * i - 2) += kk(3, 0);
    K(2 * j - 1, 2 * i - 1) += kk(3, 1);
    K(2 * j - 1, 2 * j - 2) += kk(3, 2);
    K(2 * j - 1, 2 * j - 1) += kk(3, 3);
    K(2 * j - 1, 2 * m - 2) += kk(3, 4);
    K(2 * j - 1, 2 * m - 1) += kk(3, 5);
    K(2 * m - 2, 2 * i - 2) += kk(4, 0);
    K(2 * m - 2, 2 * i - 1) += kk(4, 1);
    K(2 * m - 2, 2 * j - 2) += kk(4, 2);
    K(2 * m - 2, 2 * j - 1) += kk(4, 3);
    K(2 * m - 2, 2 * m - 2) += kk(4, 4);
    K(2 * m - 2, 2 * m - 1) += kk(4, 5);
    K(2 * m - 1, 2 * i - 2) += kk(5, 0);
    K(2 * m - 1, 2 * i - 1) += kk(5, 1);
    K(2 * m - 1, 2 * j - 2) += kk(5, 2);
    K(2 * m - 1, 2 * j - 1) += kk(5, 3);
    K(2 * m - 1, 2 * m - 2) += kk(5, 4);
    K(2 * m - 1, 2 * m - 1) += kk(5, 5);
  }
}

void Simulation::InitializeVelocity() {
  velocity = Eigen::VectorXd::Zero(mesh->rows());
}

void Simulation::InitializeAcceleration() {
  acceleration = Eigen::VectorXd::Zero(mesh->rows());
}

void Simulation::AssembleMassMatrix() {
  std::vector<Eigen::Triplet<double>> mass_entries;
  mass_entries.reserve(mesh->rows());

  M.resize(mesh->rows(), mesh->rows());

  for (int i = 0; i < mesh->rows(); ++i) {
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
