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
    float mass_, float E_, float NU_, std::shared_ptr<Mesh>& mesh_,
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

void LinearTetrahedral::Update() {}

void LinearTetrahedral::Integrate() {}

void LinearTetrahedral::AssembleForces() {
  F_ext = Eigen::VectorXf::Zero(mesh->size());
}

void LinearTetrahedral::ComputeElementStiffness(
    Eigen::Matrix12f& element_stiffness, const Eigen::Vector3f& shape_one,
    const Eigen::Vector3f& shape_two, const Eigen::Vector3f& shape_three,
    const Eigen::Vector3f& shape_four) {
  Eigen::MatrixXf B;
  AssembleStrainRelationshipMatrix(B, shape_one, shape_two, shape_three,
                                   shape_four);

  Eigen::Matrix66f D;
  AssembleStressStrainMatrix(D);

  const float V =
      ComputeElementVolume(shape_one, shape_two, shape_three, shape_four);
   element_stiffness = V * B.transpose() * D * B;
}

void LinearTetrahedral::AssembleElementStiffness() {
  for (int i = 0; i < mesh->faces_size(); i += kFaceStride) {
    std::vector<int> stiffness_coordinates;
    // Get the index face value
    int index = mesh->faces(i);
    stiffness_coordinates.push_back(index);
    const Eigen::Vector3f shape_one =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));
    index = mesh->faces(i + 1);
    stiffness_coordinates.push_back(index);
    const Eigen::Vector3f shape_two =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));

    index = mesh->faces(i + 2);
    stiffness_coordinates.push_back(index);
    const Eigen::Vector3f shape_three =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));
    index = mesh->faces(i + 3);
    stiffness_coordinates.push_back(index);
    const Eigen::Vector3f shape_four =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));

    // Eigen::Matrix12f element_stiffness_matrix;
    // ComputeElementStiffness(element_stiffness_matrix, shape_one, shape_two,
    //                        shape_three, shape_four);

    // const ElementStiffness element_stiffness = {element_stiffness_matrix,
    //                                            stiffness_coordinates};

    // k.emplace_back(element_stiffness);
  }
}

void LinearTetrahedral::AssembleElementStresses(
    Eigen::VectorXf nodal_displacement) {
  for (int i = 0; i < mesh->faces_size(); i += kFaceStride) {
    // Get the index face value
    int index = mesh->faces(i);
    const Eigen::Vector3f shape_one =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));
    index = mesh->faces(i + 1);
    const Eigen::Vector3f shape_two =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));
    index = mesh->faces(i + 2);
    const Eigen::Vector3f shape_three =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));
    index = mesh->faces(i + 3);
    const Eigen::Vector3f shape_four =
        Eigen::Vector3f(mesh->positions(index), mesh->positions(index + 1),
                        mesh->positions(index + 2));

    Eigen::MatrixXf B;
    AssembleStrainRelationshipMatrix(B, shape_one, shape_two, shape_three,
                                     shape_four);

    Eigen::Matrix66f D;
    AssembleStressStrainMatrix(D);

    const Eigen::Vector6f sigma = D * B * nodal_displacement;

    sigmas.emplace_back(sigma);
  }
}

void LinearTetrahedral::AssembleElementPlaneStresses() {
  for (const auto& sigma : sigmas) {
    const float s1 = sigma.sum();
    const float s2 =
        (sigma(0) * sigma(1) + sigma(0) * sigma(2) + sigma(1) * sigma(2)) -
        (sigma(3) * sigma(3) - sigma(4) * sigma(4) - sigma(5) * sigma(5));

    Eigen::Matrix3d ms3;
    ms3.row(0) << sigma(0), sigma(3), sigma(5);
    ms3.row(1) << sigma(3), sigma(1), sigma(4);
    ms3.row(2) << sigma(5), sigma(4), sigma(2);

    const float s3 = ms3.determinant();

    const Eigen::Vector3f plane_stress(s1, s2, s3);
    plane_stresses.emplace_back(plane_stress);
  }
}

void LinearTetrahedral::AssembleGlobalStiffness() {
  // Allocate space in K for 3nx3n elements
  const float K_rowsize = mesh->positions.size() * 3;

  K = Eigen::MatrixXf::Zero(K_rowsize, K_rowsize);

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

void LinearTetrahedral::AssembleStressStrainMatrix(Eigen::Matrix66f& D) {
  D.row(0) << 1 - NU, NU, NU, 0, 0, 0;
  D.row(1) << NU, 1 - NU, NU, 0, 0, 0;
  D.row(2) << NU, NU, 1 - NU, 0, 0, 0;
  D.row(3) << 0, 0, 0, (1 - 2 * NU) / 2, 0, 0;
  D.row(4) << 0, 0, 0, 0, (1 - 2 * NU) / 2, 0;
  D.row(5) << 0, 0, 0, 0, 0, (1 - 2 * NU) / 2;
  D *= E / ((1 + NU) * (1 - 2 * NU));
}

void LinearTetrahedral::AssembleStrainRelationshipMatrix(
    Eigen::MatrixXf& strain_relationship, const Eigen::Vector3f& shape_one,
    const Eigen::Vector3f& shape_two, const Eigen::Vector3f& shape_three,
    const Eigen::Vector3f& shape_four) {
  const auto create_beta_submatrix = [](float beta, float gamma,
                                        float delta) -> BetaSubmatrixf {
    BetaSubmatrixf B;
    B.row(0) << beta, 0, 0;
    B.row(1) << 0, gamma, 0;
    B.row(2) << 0, 0, delta;
    B.row(3) << gamma, beta, 0;
    B.row(4) << 0, delta, gamma;
    B.row(5) << delta, 0, beta;
    return B;
  };

  const float x1 = shape_one.x();
  const float y1 = shape_one.y();
  const float z1 = shape_one.z();

  const float x2 = shape_two.x();
  const float y2 = shape_two.y();
  const float z2 = shape_two.z();

  const float x3 = shape_three.x();
  const float y3 = shape_three.y();
  const float z3 = shape_three.z();

  const float x4 = shape_four.x();
  const float y4 = shape_four.y();
  const float z4 = shape_four.z();

  const float beta_1 =
      -1 * ConstructShapeFunctionParameter(y2, z2, y3, z3, y4, z4);
  const float beta_2 = ConstructShapeFunctionParameter(y1, z1, y3, z3, y4, z4);
  const float beta_3 =
      -1 * ConstructShapeFunctionParameter(y1, z1, y2, z2, y4, z4);
  const float beta_4 = ConstructShapeFunctionParameter(y1, z1, y2, z2, y3, z3);

  const float gamma_1 = ConstructShapeFunctionParameter(x2, z2, x3, z3, x4, z4);
  const float gamma_2 =
      -1 * ConstructShapeFunctionParameter(x1, z1, x3, z3, x4, z4);
  const float gamma_3 = ConstructShapeFunctionParameter(x1, z1, x2, z2, x4, z4);
  const float gamma_4 =
      -1 * ConstructShapeFunctionParameter(x1, z1, x2, z2, x3, z3);

  const float delta_1 =
      -1 * ConstructShapeFunctionParameter(x2, y2, x3, y3, x4, y4);
  const float delta_2 = ConstructShapeFunctionParameter(x1, y1, x3, y3, x4, y4);
  const float delta_3 =
      -1 * ConstructShapeFunctionParameter(x1, y1, x2, y2, x4, y4);
  const float delta_4 = ConstructShapeFunctionParameter(x1, y1, x2, y2, x3, y3);

  const BetaSubmatrixf B1 = create_beta_submatrix(beta_1, gamma_1, delta_1);
  const BetaSubmatrixf B2 = create_beta_submatrix(beta_2, gamma_2, delta_2);
  const BetaSubmatrixf B3 = create_beta_submatrix(beta_3, gamma_3, delta_3);
  const BetaSubmatrixf B4 = create_beta_submatrix(beta_4, gamma_4, delta_4);

  // Matrix is 6 x 12
  strain_relationship.resize(B1.size(), B1.cols() * 4);
  strain_relationship << B1, B2, B3, B4;
}

float LinearTetrahedral::ConstructShapeFunctionParameter(float p1, float p2,
                                                         float p3, float p4,
                                                         float p5, float p6) {
  Eigen::Matrix3d parameter;
  parameter.row(0) << 1, p1, p2;
  parameter.row(1) << 1, p3, p4;
  parameter.row(2) << 1, p5, p6;
  return parameter.determinant();
}

void LinearTetrahedral::InitializeVelocity() {
  velocity = Eigen::VectorXf::Zero(mesh->size());
}

void LinearTetrahedral::InitializeAcceleration() {
  acceleration = Eigen::VectorXf::Zero(mesh->size());
}

void LinearTetrahedral::AssembleMassMatrix() {
  // std::vector<Eigen::Triplet<float>> mass_entries;
  // mass_entries.reserve(mesh->size());

  // M.resize(mesh->size(), mesh->size());

  // for (int i = 0; i < mesh->size(); ++i) {
  //  mass_entries.push_back(Eigen::Triplet<float>(i, i, mass));
  //}

  // M.setFromTriplets(mass_entries.begin(), mass_entries.end());

  //// Set effective mass matrix
  // M_hat = (a0 * M).eval();

  //// Triangularize M_hat
  // Eigen::SimplicialLDLT<Eigen::SparseMatrixXf> solver;
  // solver.compute(M_hat);

  // assert(solver.info() == Eigen::Success && "SOLVER FAILED FOR M_HAT");

  // M_hat_inverse = solver.solve(M_hat);
}

void LinearTetrahedral::InitializeIntegrationConstants() {
  a0 = 1 / std::pow(timestep_size, 2);
  a1 = 1 / (timestep_size * 2);
  a2 = 2 * a0;
  a3 = 1 / a2;
}

Eigen::VectorXf LinearTetrahedral::SolveU(Eigen::MatrixXf k, Eigen::VectorXf f,
                                          Eigen::VectorXi indices) {
  assert(k.size() == f.size() && "K AND F DO NOT MATCH IN NUMBER OF ROWS");

  Eigen::VectorXf u;
  Eigen::VectorXf U;

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
    const float val = u(i);
    const unsigned int index = indices(i);

    U.row(index) << val;
  }

  return U;
}
/*
void LinearTetrahedral::Solve() {
  // Per-element stiffness applied to our boundary conditions
  Eigen::MatrixXf kk;

  // Our local force vector
  Eigen::VectorXf f;
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

  const Eigen::VectorXf U = SolveU(kk, f, kept_indices);

  // Set global force
  F_ext = K * U;

  // Calculate Element Stresses
  for (std::size_t i = 0; i < mesh->size(); i += kStride) {
    const float x1 = mesh->positions(i);
    const float y1 = mesh->positions(i + 1);
    const float z1 = mesh->positions(i + 2);

    const float x2 = mesh->positions(i + 3);
    const float y2 = mesh->positions(i + 4);
    const float z2 = mesh->positions(i + 5);

    const float x3 = mesh->positions(i + 6);
    const float y3 = mesh->positions(i + 7);
    const float z3 = mesh->positions(i + 8);

    const float x4 = mesh->positions(i + 9);
    const float y4 = mesh->positions(i + 10);
    const float z4 = mesh->positions(i + 11);

    // 2 * # of nodes in this segment
    Eigen::Vector12f nodal_displacement = Eigen::Vector6f::Zero();

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

float LinearTetrahedral::ComputeElementVolume(
    const Eigen::Vector3f& shape_one, const Eigen::Vector3f& shape_two,
    const Eigen::Vector3f& shape_three, const Eigen::Vector3f& shape_four) {
  float x1 = shape_one.x();
  float y1 = shape_one.y();
  float z1 = shape_one.z();

  float x2 = shape_two.x();
  float y2 = shape_two.y();
  float z2 = shape_two.z();

  float x3 = shape_three.x();
  float y3 = shape_three.y();
  float z3 = shape_three.z();

  float x4 = shape_four.x();
  float y4 = shape_four.y();
  float z4 = shape_four.z();

  Eigen::Matrix4f V;
  V.row(0) << 1, x1, y1, z1;
  V.row(1) << 1, x2, y2, z2;
  V.row(2) << 1, x3, y3, z3;
  V.row(3) << 1, x4, y4, z4;

  return V.determinant() / 6;
}
