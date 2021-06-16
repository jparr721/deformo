#include <memory>

#include "pch.h"

#define TETLIBRARY
#include "../deformo/ExplicitCentralDifference.h"
#include "../deformo/LinearTetrahedral.cpp"
#include "../deformo/LinearTetrahedral.h"
#include "../deformo/Mesh.cpp"
#include "../deformo/Mesh.h"
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
  const auto mesh = MakeBasicMesh();
  const auto lt = std::make_unique<LinearTetrahedral>(
      0.1f, 0.1f, 1.f, mesh,
      std::vector<BoundaryCondition>{{1, Eigen::Vector3f(0, 0, 0)}});

  EXPECT_TRUE(lt.get() != nullptr);
}

TEST(TestLinearTetrahedral, TestElementStiffness) {
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

  const auto lt =
      std::make_unique<LinearTetrahedral>(210e6, 0.3, 1.f, mesh, bcs);

  lt->Solve();

  Eigen::VectorXf force_compare(24);
  force_compare << -31.3293f, -5.34904f, -9.32853f, 30.7042f, -4.02576f,
      -3.07762f, -4.05312e-06f, 3.12495f, -6.70552e-07f, 1.43051e-06f, 6.24987f,
      -3.57628e-07f, -30.7042f, -4.02584f, 3.07763f, 31.3293f, -5.34919f,
      9.32853f, -3.8147e-06f, 6.24997f, -1.19209e-07f, 0.f, 3.12509f,
      9.53674e-07f;
  utils::GTestDebugPrint(lt->global_force);
  ASSERT_TRUE(lt->global_force.isApprox(force_compare));
}