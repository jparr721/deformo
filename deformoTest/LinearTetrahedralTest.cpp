#include <memory>

#include "../deformo/LinearTetrahedral.h"
#include "../deformo/Mesh.h"
#include "pch.h"

std::shared_ptr<Mesh> MakeBasicMesh() {
  Eigen::MatrixXf V(8, 3);
  V.row(0) << 0, 0, 0;
  V.row(1) << 0.025, 0, 0;
  V.row(2) << 0, 0.5, 0;
  V.row(3) << 0.025, 0.5, 0;
  V.row(4) << 0, 0, 0.25;
  V.row(5) << 0.025, 0, 0.25;
  V.row(6) << 0, 0.5, 0.25;
  V.row(7) << 0.025, 0.5, 0.25;

  Eigen::MatrixXf T(5, 4);
  T.row(0) << 1, 2, 4, 6;
  T.row(1) << 1, 4, 3, 7;
  T.row(2) << 6, 5, 7, 1;
  T.row(3) << 6, 7, 8, 4;
  T.row(4) << 1, 6, 4, 7;
}
TEST(TestLinearTetrahedral, TestConstructor) {
  const auto lt = std::make_unique<LinearTetrahedral>(0.1f, 0.1f, 1.f);
  EXPECT_TRUE(lt.get() != nullptr);
}