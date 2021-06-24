#include "pch.h"

#include <memory>

#define TETLIBRARY
#include "../deformo/ExplicitCentralDifference.h"
#include "../deformo/LinearTetrahedral.cpp"
#include "../deformo/LinearTetrahedral.h"
#include "../deformo/Mesh.cpp"
#include "../deformo/Mesh.h"
#include "../deformo/MeshGenerator.h"
#include "../deformo/MeshGenerator.cpp"
#include "tetgen.h"

auto MakeBasicMesh() -> std::shared_ptr<Mesh> {
  Eigen::MatrixXf V(8, 3);
  V.row(0) << 0, 0, 0;
  V.row(1) << 0.025, 0, 0;
  V.row(2) << 0, 0.5, 0;
  V.row(3) << 0.025, 0.5, 0;
  V.row(4) << 0, 0, 0.25;
  V.row(5) << 0.025, 0, 0.25;
  V.row(6) << 0, 0.5, 0.25;
  V.row(7) << 0.025, 0.5, 0.25;

  Eigen::MatrixXi T(5, 4);
  T.row(0) << 0, 1, 3, 5;
  T.row(1) << 0, 3, 2, 6;
  T.row(2) << 5, 4, 6, 0;
  T.row(3) << 5, 6, 7, 3;
  T.row(4) << 0, 5, 3, 6;

  return std::make_shared<Mesh>(V, T);
}

TEST(TestLinearTetrahedral, TestConstructor) {
  const float youngs_modulus = 210e6;
  const float poissons_ratio = 0.3;
  const auto mesh = MakeBasicMesh();
  const auto lt = std::make_unique<LinearTetrahedral>(
      youngs_modulus, poissons_ratio, mesh,
      std::vector<BoundaryCondition>{{1, Eigen::Vector3f(0, 0, 0)}});

  EXPECT_TRUE(lt.get() != nullptr);
}

TEST(TestLinearTetrahedral, TestElementStiffness) {
  const float youngs_modulus = 210e6;
  const float poissons_ratio = 0.3;
  const auto mesh = MakeBasicMesh();

  const auto bc_1 = BoundaryCondition{
      2 * 3,
      Eigen::Vector3f(0.f, 3.125f, 0.f),
  };

  const auto bc_2 = BoundaryCondition{
      3 * 3,
      Eigen::Vector3f(0.f, 6.25f, 0.f),
  };

  const auto bc_3 = BoundaryCondition{
      6 * 3,
      Eigen::Vector3f(0.f, 6.25f, 0.f),
  };

  const auto bc_4 = BoundaryCondition{
      7 * 3,
      Eigen::Vector3f(0.f, 3.125f, 0.f),
  };

  const std::vector bcs{{bc_1, bc_2, bc_3, bc_4}};

  const auto lt = std::make_unique<LinearTetrahedral>(
      youngs_modulus, poissons_ratio, mesh, bcs);

  const Eigen::MatrixXf plane_stresses =
      lt->Solve(youngs_modulus, poissons_ratio, mesh);

  Eigen::VectorXf p1_compare(3);
  p1_compare << 0, 0.123f, 7.4534f;
  p1_compare *= 1e9;
  ASSERT_TRUE(plane_stresses.row(0).isApprox(p1_compare));
}