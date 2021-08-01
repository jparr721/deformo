#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Homogenization.cpp"
#include "../deformo/Homogenization.h"
#include "../deformo/Numerics.h"
#include "../deformo/Rve.cpp"
#include "../deformo/Rve.h"

const auto rve = std::make_shared<Rve>(10, 10, 10, 0.1, 0.1, 0.1);

TEST(TestHomogenization, TestHexahedron) {
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(rve);
  ASSERT_TRUE(homogenization.get() != nullptr);

  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);
  ASSERT_EQ(hexahedron.size(), 4);
}

TEST(TestHomogenization, TestComputeDegreesOfFreedom) {
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(rve);
  ASSERT_TRUE(homogenization.get() != nullptr);

  const MatrixX<int> edof = homogenization->ComputeDegreesOfFreedom(1000);

  ASSERT_EQ(edof.rows(), 1000);
  ASSERT_EQ(edof.cols(), 24);
  MatrixX<int> row_0_comp(1, 24);
  row_0_comp << 4, 5, 6, 37, 38, 39, 34, 35, 36, 1, 2, 3, 367, 368, 369, 400,
      401, 402, 397, 398, 399, 364, 365, 366;
  ASSERT_TRUE(edof.row(0).isApprox(row_0_comp));

  MatrixX<int> row_last_comp(1, 24);
  row_last_comp << 3595, 3596, 3597, 3628, 3629, 3630, 3625, 3626, 3627, 3592,
      3593, 3594, 3958, 3959, 3960, 3991, 3992, 3993, 3988, 3989, 3990, 3955,
      3956, 3957;
  ASSERT_TRUE(edof.row(999).isApprox(row_last_comp));
}
