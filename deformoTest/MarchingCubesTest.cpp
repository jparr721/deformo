#include "pch.h"

#include <memory>

#include "../deformo/Numerics.h"
#include "../deformo/MarchingCubes.h"
#include "../deformo/MarchingCubes.cpp"
#include "../deformo/ImplicitSurfaceGenerator.h"

TEST(TestMarchingCubes, TestConstructor) {
  Tensor3r scalar_field;
  scalar_field.SetConstant(1);
  Real* _d = scalar_field.Instance().data();
  const auto marching_cubes = std::make_unique<MarchingCubes>(1, 1, _d);
  ASSERT_NE(marching_cubes.get(), nullptr);
}

TEST(TestMarchingCubes, TestGeneration) {
  const Material m1 = MaterialFromEandv(1, "m1", 10, 10);
  const Material m2 = MaterialFromEandv(0, "m2", 10, 10);
  const ImplicitSurfaceGenerator<Real>::Inclusion inclusion{15, 4, 5, 5, 5};
  const auto generator = std::make_unique<ImplicitSurfaceGenerator<Real>>(
      49, 49, 49,
      ImplicitSurfaceGenerator<
          Real>::ImplicitSurfaceCharacteristics::kIsotropic,
      ImplicitSurfaceGenerator<Real>::ImplicitSurfaceMicrostructure::kComposite,
      inclusion, m1.number, m2.number);

  ASSERT_NE(generator.get(), nullptr);

  Tensor3r scalar_field = generator->Generate();
  Real* _d = scalar_field.Instance().data();
  const auto marching_cubes = std::make_unique<MarchingCubes>(1, 1, _d);

  ASSERT_NE(marching_cubes.get(), nullptr);

  MatrixXr V;
  MatrixX<int> F;
  marching_cubes->GenerateGeometry(V, F, 49, 49, 49);
}
