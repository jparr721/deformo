#include "pch.h"

#include <memory>

#define TETLIBRARY
#include "../deformo/Integrators.cpp"
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
      0.1f, 0.1f, 1.f, mesh, std::vector<BoundaryCondition>{});

  EXPECT_TRUE(lt.get() != nullptr);
}

TEST(TestLinearTetrahedral, TestElementStiffness) {
  const auto mesh = MakeBasicMesh();

  const auto bc_1 = BoundaryCondition{
      2,
      Eigen::Vector3f(0.f, 3.125f, 0.f),
  };

  const auto bc_2 = BoundaryCondition{
      3,
      Eigen::Vector3f(0.f, 6.25f, 0.f),
  };

  const auto bc_3 = BoundaryCondition{
      6,
      Eigen::Vector3f(0.f, 6.25f, 0.f),
  };

  const auto bc_4 = BoundaryCondition{
      7,
      Eigen::Vector3f(0.f, 3.125f, 0.f),
  };

  const std::vector bcs{{bc_1, bc_2, bc_3, bc_4}};

  const auto lt = std::make_unique<LinearTetrahedral>(210e6, 0.3, 1.f, mesh, bcs);

  lt->Solve();
}