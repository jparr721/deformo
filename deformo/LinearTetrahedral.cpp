#define _USE_MATH_DEFINES

#include "LinearTetrahedral.h"

#include <igl/slice.h>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "Integrators.h"

LinearTetrahedral::LinearTetrahedral(
    double mass_, double E_, double NU_, std::shared_ptr<Mesh>& mesh_,
    const std::vector<BoundaryCondition>& boundary_conditions_)
    : E(E_),
      NU(NU_),
      mass(mass_),
      mesh(mesh_),
      boundary_conditions(boundary_conditions_) {
  // TODO(@jparr721) - This function should change to the integrators file once
  // we do others.
  // InitializeIntegrationConstants();

  //// Assume no initial forces so U^-dt is always the same.
  // last_displacement = mesh->positions;

  // AssembleForces();
  // AssembleElementStiffness();
  // AssembleGlobalStiffness();
  // AssembleMassMatrix();
}

void LinearTetrahedral::Update() {
}

void LinearTetrahedral::Integrate() {
}

void LinearTetrahedral::AssembleForces() {
  F_ext = Eigen::VectorXd::Zero(mesh->size());
}

void LinearTetrahedral::AssembleElementStiffness() {
  for (int i = 0; i < mesh->size(); i += kStride) {
    const double x1 = mesh->positions(i);
    const double y1 = mesh->positions(i + 1);
    const double z1 = mesh->positions(i + 2);

    const double x2 = mesh->positions(i + 3);
    const double y2 = mesh->positions(i + 4);
    const double z2 = mesh->positions(i + 5);

    const double x3 = mesh->positions(i + 6);
    const double y3 = mesh->positions(i + 7);
    const double z3 = mesh->positions(i + 8);

    const double x4 = mesh->positions(i + 9);
    const double y4 = mesh->positions(i + 10);
    const double z4 = mesh->positions(i + 11);

    const Eigen::MatrixXd B = AssembleStrainRelationshipMatrix(
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

    const Eigen::Matrix66d D = AssembleStressStrainMatrix();
    const double V =
        ComputeElementVolume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
    const Eigen::Matrix12d kk = V * B.transpose() * D * B;

    // const auto i_node_number = mesh->node_number(i / mesh->kNumDimensions);
    // const auto j_node_number =
    //    mesh->node_number((i + 3) / mesh->kNumDimensions);
    // const auto m_node_number =
    //    mesh->node_number((i + 6) / mesh->kNumDimensions);
    // const auto n_node_number =
    //    mesh->node_number((i + 9) / mesh->kNumDimensions);

    // const ElementStiffness element_stiffness = {
    //    kk, std::vector<unsigned int>{i_node_number, j_node_number,
    //                                  m_node_number, n_node_number}};

    //k.emplace_back(element_stiffness);
  }
}

/*
void LinearTetrahedral::AssembleElementStresses(Eigen::VectorXd
nodal_displacement) { for (int i = 0; i < mesh->size(); i += kStride) { const
double x1 = mesh->positions(i); const double y1 = mesh->positions(i + 1); const
double z1 = mesh->positions(i + 2);

    const double x2 = mesh->positions(i + 3);
    const double y2 = mesh->positions(i + 4);
    const double z2 = mesh->positions(i + 5);

    const double x3 = mesh->positions(i + 6);
    const double y3 = mesh->positions(i + 7);
    const double z3 = mesh->positions(i + 8);

    const double x4 = mesh->positions(i + 9);
    const double y4 = mesh->positions(i + 10);
    const double z4 = mesh->positions(i + 11);

    const Eigen::MatrixXd B = AssembleStrainRelationshipMatrix(
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

    const Eigen::Matrix66d D = AssembleStressStrainMatrix();

    const Eigen::Vector6d sigma = D * B * nodal_displacement;

    sigmas.emplace_back(sigma);
  }
}

void LinearTetrahedral::AssembleElementPlaneStresses() {
  for (const auto& sigma : sigmas) {
    const double s1 = sigma.sum();
    const double s2 =
        (sigma(0) * sigma(1) + sigma(0) * sigma(2) + sigma(1) * sigma(2)) -
        (sigma(3) * sigma(3) - sigma(4) * sigma(4) - sigma(5) * sigma(5));

    Eigen::Matrix3d ms3;
    ms3.row(0) << sigma(0), sigma(3), sigma(5);
    ms3.row(1) << sigma(3), sigma(1), sigma(4);
    ms3.row(2) << sigma(5), sigma(4), sigma(2);

    const double s3 = ms3.determinant();

    const Eigen::Vector3d plane_stress(s1, s2, s3);
    plane_stresses.emplace_back(plane_stress);
  }
}

void LinearTetrahedral::AssembleGlobalStiffness() {
  // Allocate space in K for 3nx3n elements
  const double K_rowsize = mesh->positions.size() * mesh->kNumDimensions;

  K = Eigen::MatrixXd::Zero(K_rowsize, K_rowsize);

  // Indices in the element stiffness matrix to collect
  int kk_0 = -1;
  int kk_1 = 0;

  // Iterate through all points and element stiffness matrices, generating the
  // output inside of the K matrix. Lotta, shallow loops, O(144)
  for (const auto& element_stiffness : k) {
    const auto kk = element_stiffness.stiffness_matrix;
    const auto pts = element_stiffness.indices;

    for (const auto& pt : pts) {
      for (int i = 3; i > 0; --i) {
        const int left_K_idx = 3 * pt - i;
        kk_0 += 1;

        for (const auto& ppt : pts) {
          for (int j = 3; j > 0; --j) {
            // If our right-submatrix coordinate (y) > the cols of the
            // submatrix, reset
            if (kk_1 == kk.cols()) {
              kk_1 = 0;
            }

            const int right_K_idx = 3 * ppt - j;

            K(left_K_idx, right_K_idx) += kk(kk_0, kk_1);
          }
        }
      }
    }
  }
}
*/
Eigen::Matrix66d LinearTetrahedral::AssembleStressStrainMatrix() {
  Eigen::Matrix66d D;
  D.row(0) << 1 - NU, NU, NU, 0, 0, 0;
  D.row(1) << NU, 1 - NU, NU, 0, 0, 0;
  D.row(2) << NU, NU, 1 - NU, 0, 0, 0;
  D.row(3) << 0, 0, 0, (1 - 2 * NU) / 2, 0, 0;
  D.row(4) << 0, 0, 0, 0, (1 - 2 * NU) / 2, 0;
  D.row(5) << 0, 0, 0, 0, 0, (1 - 2 * NU) / 2;

  D *= E / ((1 + NU) * (1 - 2 * NU));
  return D;
}

Eigen::MatrixXd LinearTetrahedral::AssembleStrainRelationshipMatrix(
    double x1, double y1, double z1, double x2, double y2, double z2, double x3,
    double y3, double z3, double x4, double y4, double z4) {
  const double beta_1 =
      -1 * ConstructShapeFunctionParameter(y2, z2, y3, z3, y4, z4);
  const double beta_2 = ConstructShapeFunctionParameter(y1, z1, y3, z3, y4, z4);
  const double beta_3 =
      -1 * ConstructShapeFunctionParameter(y1, z1, y2, z2, y4, z4);
  const double beta_4 = ConstructShapeFunctionParameter(y1, z1, y2, z2, y3, z3);

  const double gamma_1 =
      ConstructShapeFunctionParameter(x2, z2, x3, z3, x4, z4);
  const double gamma_2 =
      -1 * ConstructShapeFunctionParameter(x1, z1, x3, z3, x4, z4);
  const double gamma_3 =
      ConstructShapeFunctionParameter(x1, z1, x2, z2, x4, z4);
  const double gamma_4 =
      -1 * ConstructShapeFunctionParameter(x1, z1, x2, z2, x3, z3);

  const double delta_1 =
      -1 * ConstructShapeFunctionParameter(x2, y2, x3, y3, x4, y4);
  const double delta_2 =
      ConstructShapeFunctionParameter(x1, y1, x3, y3, x4, y4);
  const double delta_3 =
      -1 * ConstructShapeFunctionParameter(x1, y1, x2, y2, x4, y4);
  const double delta_4 =
      ConstructShapeFunctionParameter(x1, y1, x2, y2, x3, y3);

  BetaSubmatrixd B1;
  B1.row(0) << beta_1, 0, 0;
  B1.row(1) << 0, gamma_1, 0;
  B1.row(2) << 0, 0, delta_1;
  B1.row(3) << gamma_1, beta_1, 0;
  B1.row(4) << 0, delta_1, gamma_1;
  B1.row(5) << delta_1, 0, beta_1;

  BetaSubmatrixd B2;
  B2.row(0) << beta_2, 0, 0;
  B2.row(1) << 0, gamma_2, 0;
  B2.row(2) << 0, 0, delta_2;
  B2.row(3) << gamma_2, beta_2, 0;
  B2.row(4) << 0, delta_2, gamma_2;
  B2.row(5) << delta_2, 0, beta_2;

  BetaSubmatrixd B3;
  B3.row(0) << beta_3, 0, 0;
  B3.row(1) << 0, gamma_3, 0;
  B3.row(2) << 0, 0, delta_3;
  B3.row(3) << gamma_3, beta_3, 0;
  B3.row(4) << 0, delta_3, gamma_3;
  B3.row(5) << delta_3, 0, beta_3;

  BetaSubmatrixd B4;
  B4.row(0) << beta_4, 0, 0;
  B4.row(1) << 0, gamma_4, 0;
  B4.row(2) << 0, 0, delta_4;
  B4.row(3) << gamma_4, beta_4, 0;
  B4.row(4) << 0, delta_4, gamma_4;
  B4.row(5) << delta_4, 0, beta_4;

  // Matrix is 6 x 12
  Eigen::MatrixXd B(B1.size(), B1.cols() * 4);
  B << B1, B2, B3, B4;

  return B;
}

double LinearTetrahedral::ConstructShapeFunctionParameter(double p1, double p2,
                                                          double p3, double p4,
                                                          double p5,
                                                          double p6) {
  Eigen::Matrix3d parameter;
  parameter.row(0) << 1, p1, p2;
  parameter.row(1) << 1, p3, p4;
  parameter.row(2) << 1, p5, p6;
  return parameter.determinant();
}

void LinearTetrahedral::InitializeVelocity() {
  velocity = Eigen::VectorXd::Zero(mesh->size());
}

void LinearTetrahedral::InitializeAcceleration() {
  acceleration = Eigen::VectorXd::Zero(mesh->size());
}

void LinearTetrahedral::AssembleMassMatrix() {
  //std::vector<Eigen::Triplet<double>> mass_entries;
  //mass_entries.reserve(mesh->size());

  //M.resize(mesh->size(), mesh->size());

  //for (int i = 0; i < mesh->size(); ++i) {
  //  mass_entries.push_back(Eigen::Triplet<double>(i, i, mass));
  //}

  //M.setFromTriplets(mass_entries.begin(), mass_entries.end());

  //// Set effective mass matrix
  //M_hat = (a0 * M).eval();

  //// Triangularize M_hat
  //Eigen::SimplicialLDLT<Eigen::SparseMatrixXd> solver;
  //solver.compute(M_hat);

  //assert(solver.info() == Eigen::Success && "SOLVER FAILED FOR M_HAT");

  //M_hat_inverse = solver.solve(M_hat);
}

void LinearTetrahedral::InitializeIntegrationConstants() {
  a0 = 1 / std::pow(timestep_size, 2);
  a1 = 1 / (timestep_size * 2);
  a2 = 2 * a0;
  a3 = 1 / a2;
}

Eigen::VectorXd LinearTetrahedral::SolveU(Eigen::MatrixXd k, Eigen::VectorXd f,
                                          Eigen::VectorXi indices) {
  assert(k.size() == f.size() && "K AND F DO NOT MATCH IN NUMBER OF ROWS");

  Eigen::VectorXd u;
  Eigen::VectorXd U;

  // Indices is always # of nodes not including multi-coordinate layouts
  u.resize(mesh->positions.size());

  // Begin the process of initializing global displacement
  U.setZero(u.size());

  // Full pivot LU factorization of k minimized as far as it can go,
  // then we solve with respect to f, assigning to our global displacement
  // vector.
  u = k.fullPivLu().solve(f);

  assert(indices.size() == u.size() && "INDEX AND U DIFFER");

  for (int i = 0; i < u.size(); ++i) {
    const double val = u(i);
    const unsigned int index = indices(i);

    U.row(index) << val;
  }

  return U;
}
/*
void LinearTetrahedral::Solve() {
  // Per-element stiffness applied to our boundary conditions
  Eigen::MatrixXd kk;

  // Our local force vector
  Eigen::VectorXd f;
  // Resize the force vector to the # of boundary conditions * 2;
  f.resize(boundary_conditions.size() * 2);

  // Stacked vector of xy columns/rows to slice.
  Eigen::VectorXi kept_indices;
  kept_indices.resize(boundary_conditions.size() * 2);

  // So we have a key difference from literature. Our boundary conditions
  // specify external forces as opposed to specifying which nodes are fixed, we
  // just assume nodes without forces are fixed at 0. So in a problem of the
  // form U1 - U4 where we have 4 nodes, and 2 (U2, U3) are degrees of freedom,
  // we can make conditions as:
  //
  // { node: 2, force(x, y, z), node: 3, force(x, y, z) }
  //
  // Where x, y and z are the directions in which the force is applied. We can
map
  // into the force vector at these coordinates
  int segment = 0;  // Keep track of kept indices as a 0-indexed vector
  for (const auto& boundary_condition : boundary_conditions) {
    // Get the node number so we can begin indexing
    const unsigned int node_number = boundary_condition.node;

    // The row is the same as the index segment
    f.row(segment) << boundary_condition.force.x;
    f.row(segment + 1) << boundary_condition.force.y;
    f.row(segment + 2) << boundary_condition.force.z;

    const unsigned int k_col_x = node_number * 3 - 3;
    const unsigned int k_col_y = node_number * 3 - 2;
    const unsigned int k_col_z = node_number * 3 - 1;
    kept_indices.segment(segment, 3) << k_col_x, k_col_y, k_col_z;
    segment += 3;
  }

  igl::slice(K, kept_indices, kept_indices, kk);

  const Eigen::VectorXd U = SolveU(kk, f, kept_indices);

  // Set global force
  F_ext = K * U;

  // Calculate Element Stresses
  for (std::size_t i = 0; i < mesh->size(); i += kStride) {
    const double x1 = mesh->positions(i);
    const double y1 = mesh->positions(i + 1);
    const double z1 = mesh->positions(i + 2);

    const double x2 = mesh->positions(i + 3);
    const double y2 = mesh->positions(i + 4);
    const double z2 = mesh->positions(i + 5);

    const double x3 = mesh->positions(i + 6);
    const double y3 = mesh->positions(i + 7);
    const double z3 = mesh->positions(i + 8);

    const double x4 = mesh->positions(i + 9);
    const double y4 = mesh->positions(i + 10);
    const double z4 = mesh->positions(i + 11);

    // 2 * # of nodes in this segment
    Eigen::Vector12d nodal_displacement = Eigen::Vector6d::Zero();

    const auto i_node_number = mesh->node_number(i / mesh->kNumDimensions);
    const auto j_node_number =
        mesh->node_number((i + 3) / mesh->kNumDimensions);
    const auto m_node_number =
        mesh->node_number((i + 6) / mesh->kNumDimensions);
    const auto n_node_number =
        mesh->node_number((i + 9) / mesh->kNumDimensions);

    // Corresponds to node index in the U vector
    const unsigned int i_row_l = i_node_number * 3 - 3;
    const unsigned int i_row_m = i_node_number * 3 - 2;
    const unsigned int i_row_r = i_node_number * 3 - 1;

    const unsigned int j_row_l = j_node_number * 3 - 3;
    const unsigned int j_row_m = j_node_number * 3 - 2;
    const unsigned int j_row_r = j_node_number * 3 - 1;

    const unsigned int m_row_l = m_node_number * 3 - 3;
    const unsigned int m_row_m = m_node_number * 3 - 2;
    const unsigned int m_row_r = m_node_number * 3 - 1;

    const unsigned int n_row_l = n_node_number * 3 - 3;
    const unsigned int n_row_m = n_node_number * 3 - 2;
    const unsigned int n_row_r = n_node_number * 3 - 1;

    nodal_displacement.segment(0, 3) << i_row_l, i_row_m, i_row_r;
    nodal_displacement.segment(3, 3) << j_row_l, j_row_m, j_row_r;
    nodal_displacement.segment(6, 3) << m_row_l, m_row_m, m_row_r;
    nodal_displacement.segment(9, 3) << n_row_l, n_row_m, n_row_r;

    AssembleElementStresses(nodal_displacement);
  }
}
*/

double LinearTetrahedral::ComputeElementVolume(double x1, double y1, double z1,
                                               double x2, double y2, double z2,
                                               double x3, double y3, double z3,
                                               double x4, double y4,
                                               double z4) {
  Eigen::Matrix4d V;
  V.row(0) << 1, x1, y1, z1;
  V.row(1) << 1, x2, y2, z2;
  V.row(2) << 1, x3, y3, z3;
  V.row(3) << 1, x4, y4, z4;

  return V.determinant() / 6;
}
