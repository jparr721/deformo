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

TEST(TestNumerics, TestHStack) {
  MatrixX<int> one(2, 2);
  one.setConstant(1);

  MatrixX<int> two(2, 2);
  two.setConstant(2);

  MatrixX<int> three(2, 2);
  three.setConstant(3);

  MatrixX<int> comp(6, 2);
  comp << 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3;

  const MatrixX<int> stacked =
      linear_algebra::HStack(std::vector{one, two, three});

  ASSERT_TRUE(stacked.isApprox(comp));
}

TEST(TestNumerics, TestHStackVectors) {
  VectorX<int> one(2);
  one.setConstant(1);

  VectorX<int> two(2);
  two.setConstant(2);

  VectorX<int> three(2);
  three.setConstant(3);

  MatrixX<int> comp(3, 2);
  comp << 1, 1, 2, 2, 3, 3;

  const MatrixX<int> stacked =
      linear_algebra::HStack(std::vector{one, two, three});

  ASSERT_TRUE(stacked.isApprox(comp));
}

TEST(TestNumerics, TestSetLayer) { 
  Tensor3i t(2, 2, 1);
  t.SetConstant(1);
  MatrixX<int> layer(2, 2);
  layer.setConstant(5);

  ASSERT_TRUE(t(0, 0, 0) = 1);
  ASSERT_TRUE(t(1, 0, 0) = 1);
  ASSERT_TRUE(t(0, 1, 0) = 1);
  ASSERT_TRUE(t(1, 1, 0) = 1);

  t.SetLayer(0, layer);
  ASSERT_TRUE(t(0, 0, 0) = 5);
  ASSERT_TRUE(t(1, 0, 0) = 5);
  ASSERT_TRUE(t(0, 1, 0) = 5);
  ASSERT_TRUE(t(1, 1, 0) = 5);
}

TEST(TestNumerics, TestReArrange) { 
  MatrixX<int> m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  VectorX<int> indices(3);
  indices << 1, 2, 0;

  const MatrixX<int> m2 = linear_algebra::ReArrange(m, indices);

  MatrixX<int> comp(3, 3);
  comp.row(0) << 4, 5, 6;
  comp.row(1) << 7, 8, 9;
  comp.row(2) << 1, 2, 3;

  ASSERT_TRUE(m2.isApprox(comp));
}

TEST(TestNumerics, TestReArrange2) { 
  MatrixX<int> m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  VectorX<int> indices(4);
  indices << 1, 2, 0, 3;

  const MatrixX<int> m2 = linear_algebra::ReArrange(m, indices);

  MatrixX<int> comp(4, 3);
  comp.row(0) << 4, 5, 6;
  comp.row(1) << 7, 8, 9;
  comp.row(2) << 1, 2, 3;
  comp.row(3) << 10, 11, 12;

  ASSERT_TRUE(m2.isApprox(comp));
}

TEST(TestNumerics, TestIndexMatrixByMatrix) {
  MatrixX<int> m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  MatrixX<int> indices(3, 3);
  indices.row(0) << 0, 1, 2;
  indices.row(1) << 1, 2, 0;
  indices.row(2) << 4, 5, 2;

  const MatrixX<int> m2 = linear_algebra::IndexMatrixByMatrix(m, indices);

  MatrixX<int> comp(9, 3);
  comp.row(0) << 1, 2, 3;
  comp.row(1) << 4, 5, 6;
  comp.row(2) << 7, 8, 9;

  comp.row(3) << 4, 5, 6;
  comp.row(4) << 7, 8, 9;
  comp.row(5) << 1, 2, 3;

  comp.row(6) << 13, 14, 15;
  comp.row(7) << 16, 17, 18;
  comp.row(8) << 7, 8, 9;

  ASSERT_TRUE(m2.isApprox(comp));
}

TEST(TestNumerics, TestIndexMatrixByMatrix2) {
  MatrixXr m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  MatrixX<int> indices(3, 3);
  indices.row(0) << 0, 1, 2;
  indices.row(1) << 1, 2, 0;
  indices.row(2) << 4, 5, 2;

  const MatrixXr m2 = linear_algebra::IndexMatrixByMatrix(m, indices);

  MatrixXr comp(9, 3);
  comp.row(0) << 1, 2, 3;
  comp.row(1) << 4, 5, 6;
  comp.row(2) << 7, 8, 9;

  comp.row(3) << 4, 5, 6;
  comp.row(4) << 7, 8, 9;
  comp.row(5) << 1, 2, 3;

  comp.row(6) << 13, 14, 15;
  comp.row(7) << 16, 17, 18;
  comp.row(8) << 7, 8, 9;

  ASSERT_TRUE(m2.isApprox(comp));
}

TEST(TestNumerics, TestIndexMatrixByMatrixWithCol) {
  MatrixX<int> m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  MatrixX<int> indices(3, 3);
  indices.row(0) << 0, 1, 2;
  indices.row(1) << 1, 2, 0;
  indices.row(2) << 4, 5, 2;

  const MatrixX<int> m2 = linear_algebra::IndexMatrixByMatrix(m, indices, 0);

  MatrixX<int> comp(3, 3);
  comp.row(0) << 1, 4, 7;
  comp.row(1) << 4, 7, 1;
  comp.row(2) << 13, 16, 7;

  ASSERT_TRUE(m2.isApprox(comp));
}

TEST(TestNumerics, TestIndexMatrixByMatrixWithCol2) {
  MatrixXr m(6, 3); 
  m.row(0) << 1, 2, 3;
  m.row(1) << 4, 5, 6;
  m.row(2) << 7, 8, 9;
  m.row(3) << 10, 11, 12;
  m.row(4) << 13, 14, 15;
  m.row(5) << 16, 17, 18;

  MatrixX<int> indices(3, 3);
  indices.row(0) << 0, 1, 2;
  indices.row(1) << 1, 2, 0;
  indices.row(2) << 4, 5, 2;

  const MatrixXr m2 = linear_algebra::IndexMatrixByMatrix(m, indices, 0);

  MatrixXr comp(3, 3);
  comp.row(0) << 1, 4, 7;
  comp.row(1) << 4, 7, 1;
  comp.row(2) << 13, 16, 7;

  ASSERT_TRUE(m2.isApprox(comp));
}

