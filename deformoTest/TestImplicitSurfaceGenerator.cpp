#include "pch.h"

#include <memory>

#include "../deformo/Numerics.h"
#include "../deformo/Numerics.cpp"
#include "../deformo/ImplicitSurfaceGenerator.h"
#include "../deformo/ImplicitSurfaceGenerator.cpp"

TEST(TestImplicitSurfaceGenerator, TestConstructor) {
  const auto generator = std::make_unique<ImplicitSurfaceGenerator>(10, 10, 10);
  GTEST_ASSERT_NE(generator.get(), nullptr);
}

TEST(TestImplicitSurfaceGenerator, TestGenerator) {
  const auto generator = std::make_unique<ImplicitSurfaceGenerator>(50, 50, 50);
  const BinaryInclusion inclusion{15, 4, 5, 5, 5};

  const Tensor3r generated = generator->Generate(inclusion);
  std::cout << generated.At(0) << std::endl;
}