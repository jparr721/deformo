#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/EigenTypes.h"
#include "../deformo/ExplicitCentralDifference.cpp"
#include "../deformo/ExplicitCentralDifference.h"
#include "../deformo/LinearTetrahedral.h"
#include "../deformo/Rayleigh.cpp"
#include "../deformo/Rayleigh.h"

auto MakeStiffnessMatrix() -> Eigen::Matrix2f {
  Eigen::Matrix2f stiffness;
  stiffness << 6, -2, -2, 4;

  return stiffness;
}

auto MakeMassMatrix() -> SparseMatrixXr {
  using T = Eigen::Triplet<Real>;
  SparseMatrixXr mass_matrix;
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
  const Eigen::Vector2f initial_displacement = Eigen::Vector2f(0, 0);
  const Eigen::Vector2f initial_forces = Eigen::Vector2f(0.f, 10.f);
  const Eigen::Matrix2f stiffness = MakeStiffnessMatrix();
  const SparseMatrixXr mass_matrix = MakeMassMatrix();

  const auto integrator = std::make_unique<ExplicitCentralDifferenceMethod>(
      .28, 1, stiffness, initial_displacement, initial_forces);

  std::cout << integrator->Acceleration() << std::endl;
  ASSERT_TRUE(integrator->Acceleration().isApprox(Eigen::Vector2f(0, 10)));

  GTEST_ASSERT_NE(integrator.get(), nullptr);
}

TEST(TestExplicitCentralDifference, TestSolver) {
  VectorXr displacement = Eigen::Vector2f(0.f, 0.f);
  Eigen::Vector2f forces = Eigen::Vector2f(0.f, 10.f);
  const Eigen::Matrix2f stiffness = MakeStiffnessMatrix();
  const SparseMatrixXr mass_matrix = MakeMassMatrix();

  const auto integrator = std::make_unique<ExplicitCentralDifferenceMethod>(
      .28, mass_matrix, stiffness, displacement, forces);
  integrator->SetDamping(0, 0);

  for (int i = 0; i < 12; ++i) {
    integrator->Solve(displacement, forces);
    std::cerr << "displacement: " << std::endl;
    std::cerr << displacement << std::endl;
    std::cerr << "===" << std::endl;
    if (i == 0) {
      VectorXr compare(2);
      compare << 0, 0.392f;
      ASSERT_TRUE(compare.isApprox(displacement));
    }
  }

  VectorXr compare(2);
  compare << 1.0223f, 2.60083f;
  ASSERT_TRUE(compare.isApprox(displacement));
}
