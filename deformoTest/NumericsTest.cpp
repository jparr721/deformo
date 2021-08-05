#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Numerics.cpp"
#include "../deformo/Numerics.h"

TEST(TestNumerics, TestTensor3Constructor) {
  const Tensor3i t(3, 2, 3);
  ASSERT_TRUE(t.Dimensions().size() == 3);
  ASSERT_TRUE(t.Dimension(0) == 3);
  ASSERT_TRUE(t.Dimension(1) == 2);
  ASSERT_TRUE(t.Dimension(2) == 3);
}

TEST(TestNumerics, TestTensor3Operator) {
  Tensor3i t(2, 2, 1);
  t.SetConstant(0);
  t(0, 0, 0) = 1;
  ASSERT_TRUE(t(0, 0, 0) == 1);
}

TEST(TestNumerics, TestTensor3AppendCol) {
  Tensor3i t(2, 2, 1);
  t.SetConstant(1);

  VectorX<int> c1(2);
  c1 << 2, 3;

  std::vector<VectorX<int>> cols{c1};

  Tensor3i t2 = t.Append(cols, Tensor3i::InsertOpIndex::kStart,
                         Tensor3i::OpOrientation::kCol);
  ASSERT_TRUE(t2.Dimension(0) == 2);
  ASSERT_TRUE(t2.Dimension(1) == 3);
  ASSERT_TRUE(t2.Dimension(2) == 1);

  MatrixX<int> layer_comp(2, 3);
  layer_comp.row(0) << 2, 1, 1;
  layer_comp.row(1) << 3, 1, 1;

  const MatrixX<int> m = t2.Matrix(2, 3);

  ASSERT_TRUE(m.isApprox(layer_comp));
}

TEST(TestNumerics, TestTensor3AppendRow) {
  Tensor3i t(2, 2, 1);
  t.SetConstant(1);

  VectorX<int> c1(2);
  c1 << 2, 3;

  std::vector<VectorX<int>> cols{c1};

  Tensor3i t2 = t.Append(cols, Tensor3i::InsertOpIndex::kStart,
                         Tensor3i::OpOrientation::kRow);
  ASSERT_TRUE(t2.Dimension(0) == 3);
  ASSERT_TRUE(t2.Dimension(1) == 2);
  ASSERT_TRUE(t2.Dimension(2) == 1);

  MatrixX<int> layer_comp(3, 2);
  layer_comp.row(0) << 2, 3;
  layer_comp.row(1) << 1, 1;
  layer_comp.row(2) << 1, 1;

  const MatrixX<int> m = t2.Matrix(3, 2);

  ASSERT_TRUE(m.isApprox(layer_comp));
}

TEST(TestNumerics, TestTensor3AppendLayer) {
  Tensor3i t(2, 2, 1);
  t.SetConstant(1);

  MatrixX<int> layer(2, 2);
  layer.setConstant(2);

  Tensor3i t2 = t.Append(layer, Tensor3i::InsertOpIndex::kEnd);
  ASSERT_TRUE(t2.Dimension(0) == 2);
  ASSERT_TRUE(t2.Dimension(1) == 2);
  ASSERT_TRUE(t2.Dimension(2) == 2);

  MatrixX<int> l0_comp(2, 2);
  l0_comp.setConstant(1);
  MatrixX<int> layer_0 = t2.At(0);

  ASSERT_TRUE(layer_0.isApprox(l0_comp));

  MatrixX<int> l1_comp(2, 2);
  l1_comp.setConstant(2);
  MatrixX<int> layer_1 = t2.At(1);
  ASSERT_TRUE(layer_1.isApprox(l1_comp));
}

TEST(TestNumerics, TestIntoVector) {
  Tensor3i t(2, 2, 1);
  t.SetConstant(2);

  const VectorX<int> v = t.Vector(4);
  ASSERT_TRUE(v.rows() == 4);
}

TEST(TestNumerics, TestVStack) {
  MatrixX<int> one(2, 2);
  one.setConstant(1);

  MatrixX<int> two(2, 2);
  two.setConstant(2);

  MatrixX<int> three(2, 2);
  three.setConstant(3);

  MatrixX<int> comp(6, 2);
  comp << 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3;

  const MatrixX<int> stacked =
      linear_algebra::VStack(std::vector{one, two, three});

  ASSERT_TRUE(stacked.isApprox(comp));
}
