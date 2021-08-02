#include "pch.h"

#include <memory>

#define TETLIBRARY
#include "../deformo/ExplicitCentralDifference.h"
#include "../deformo/LinearTetrahedral.cpp"
#include "../deformo/LinearTetrahedral.h"
#include "../deformo/Mesh.cpp"
#include "../deformo/Mesh.h"
#include "../deformo/MeshGenerator.cpp"
#include "../deformo/MeshGenerator.h"
#include "tetgen.h"

auto MakeBasicMesh() -> std::shared_ptr<Mesh> {
  MatrixXr V(8, 3);
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
  const Real youngs_modulus = 210e6;
  const Real poissons_ratio = 0.3;
  const auto mesh = MakeBasicMesh();
  const auto lt = std::make_unique<LinearTetrahedral>(
      youngs_modulus, poissons_ratio, mesh,
      std::vector<BoundaryCondition>{{1, Eigen::Vector3f(0, 0, 0)}});

  EXPECT_TRUE(lt.get() != nullptr);
}

TEST(TestLinearTetrahedral, TestElementStiffness) {
  const Real youngs_modulus = 210e6;
  const Real poissons_ratio = 0.3;
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

  const MatrixXr stresses =
      lt->SolveStatic(youngs_modulus, poissons_ratio, mesh);

  MatrixXr stress_compare;
  stress_compare.resize(5, 6);
  stress_compare.row(0) << 1.47278, 3.43648, 1.47278, -0.0205161, 0.00896624, 0;
  stress_compare.row(1) << 0.00639241, 2.7694, 0.710224, -0.0128633, 0.0133638,
      -0.0703542;
  stress_compare.row(2) << 1.47278, 3.43648, 1.47278, 0.0205193, -.00896781, 0;
  stress_compare.row(3) << 0.00639215, 2.7694, 0.710224, 0.0128691, -.0133625,
      -0.0703543;
  stress_compare.row(4) << 0.00963581, 2.79405, 0.794476, -4.88281e-06,
      -1.38283e-07, 0.220376;

  ASSERT_TRUE(stresses.isApprox(stress_compare * 1e3));
}