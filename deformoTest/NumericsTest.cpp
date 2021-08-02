#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Numerics.cpp"
#include "../deformo/Numerics.h"

TEST(TestNumerics, TestTensor3) {
  const Tensor3i t(3, 3, 3);
  ASSERT_TRUE(t.Dimensions().size() == 3);
}
