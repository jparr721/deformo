#include "pch.h"
#include <memory>
#include <iostream>
#include <Eigen/SparseLU>

#include "../deformo/EigenTypes.h"
#include "../deformo/Integrators.h"
#include "../deformo/LinearTetrahedral.h"
auto MakeInitialDisplacement() -> Eigen::Vector2f {
  return Eigen::Vector2f(0.f, 0.392f);
}

auto MakeStiffnessMatrix() -> Eigen::Matrix2f {
  Eigen::Matrix2f stiffness;
  stiffness << 6, -2, -2, 4;

  return stiffness;
}

auto MakeMassMatrix() -> Eigen::SparseMatrixXf {
  using T = Eigen::Triplet<float>;
  Eigen::SparseMatrixXf mass_matrix;
  mass_matrix.resize(2, 2);
  const auto ul = T(0, 0, 2);
  const auto ur = T(0, 1, 0);
  const auto ll = T(1, 0, 0);
  const auto lr = T(1, 1, 1);
  auto vals = std::vector{ul, ur, ll, lr};
  mass_matrix.setFromTriplets(vals.begin(), vals.end());
  return mass_matrix;
}

TEST(TestExplicitCentralDifference, TestConstructor) {
  const Eigen::Vector2f initial_displacement = MakeInitialDisplacement();
  const Eigen::Matrix2f stiffness = MakeStiffnessMatrix();
  const Eigen::SparseMatrixXf mass_matrix = MakeMassMatrix();

  const auto integrator = std::make_unique<ExplicitCentralDifferenceMethod>(
      0.28, initial_displacement, stiffness, mass_matrix);

  GTEST_ASSERT_NE(integrator.get(), nullptr);
}

TEST(TestExplicitCentralDifference, TestSolver) {
  const Eigen::Vector2f initial_displacement = MakeInitialDisplacement();
  const Eigen::Matrix2f stiffness = MakeStiffnessMatrix();
  const Eigen::SparseMatrixXf mass_matrix = MakeMassMatrix();

  const auto integrator = std::make_unique<ExplicitCentralDifferenceMethod>(
      .28, initial_displacement, stiffness, mass_matrix);

  Eigen::VectorXf positions(2);
  positions << 0, 0;
  const Eigen::Vector2f force(0, 10);

  for (int i = 0; i < 12; ++i) {
    integrator->Solve(positions, force);
    if (i == 0) {
      Eigen::VectorXf compare(2);
      compare << 0, 0.392f;
     ASSERT_TRUE(compare.isApprox(positions));
    }
  }

  Eigen::VectorXf compare(2);
  compare << 1.0223f, 2.60083f;
 ASSERT_TRUE(compare.isApprox(positions));
}
