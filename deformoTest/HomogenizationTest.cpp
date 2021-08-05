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

  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(1000);

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

TEST(TestHomogenizations, TestComputeUniqueNodes) {
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(rve);
  ASSERT_TRUE(homogenization.get() != nullptr);

  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(1000);

  const MatrixX<int> first_layer = unique_nodes.At(0);
  const MatrixX<int> last_layer = unique_nodes.At(0);

  ASSERT_TRUE(first_layer.isApprox(last_layer));

  const int rows = unique_nodes.Dimension(0);
  const int cols = unique_nodes.Dimension(1);
  const int layers = unique_nodes.Dimension(2);

  // Ensure the mirroring worked as intended.
  for (int l = 0; l < layers; ++l) {
    const VectorX<int> top_row = unique_nodes.Row(l, 0);
    const VectorX<int> bottom_row = unique_nodes.Row(l, rows - 1);

    ASSERT_TRUE(top_row.isApprox(bottom_row));

    const VectorX<int> left_col = unique_nodes.Col(l, 0);
    const VectorX<int> right_col = unique_nodes.Col(l, cols - 1);

    ASSERT_TRUE(left_col.isApprox(right_col));
  }
}

TEST(TestHomogenization, TestComputeUniqueDegreesOfFreedom) {
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(rve);
  ASSERT_TRUE(homogenization.get() != nullptr);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);

  auto row_0_comp = VectorX<int>(24);
  row_0_comp << 4, 5, 6, 34, 35, 36, 31, 32, 33, 1, 2, 3, 304, 305, 306, 334,
      335, 336, 331, 332, 333, 301, 302, 303;

  ASSERT_TRUE(row_0_comp.transpose().isApprox(unique_dof.row(0)));

  auto row_n_comp = VectorX<int>(24);
  row_n_comp << 2971, 2972, 2973, 2701, 2702, 2703, 2728, 2729, 2730, 2998,
      2999, 3000, 271, 272, 273, 1, 2, 3, 28, 29, 30, 298, 299, 300;

  ASSERT_TRUE(
      row_n_comp.transpose().isApprox(unique_dof.row(unique_dof.rows() - 1)));
}
