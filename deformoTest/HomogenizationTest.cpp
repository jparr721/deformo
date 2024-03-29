#include "pch.h"

#include <Eigen/SparseLU>
#include <iostream>
#include <memory>

#include "../deformo/Homogenization.cpp"
#include "../deformo/Homogenization.h"
#include "../deformo/Numerics.h"
#include "../deformo/Rve.cpp"
#include "../deformo/Rve.h"

auto IsApprox(Real lhs, Real rhs, Real epsilon) -> bool {
  return std::fabs(lhs - rhs) < epsilon;
}

TEST(TestHomogenization, TestHexahedron) {
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), Material{1},
                                   Material{0});

  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);

  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);
  ASSERT_EQ(hexahedron.size(), 4);
}

TEST(TestHomogenization, TestComputeDegreesOfFreedom) {
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), Material{1},
                                   Material{0});
  rve->ComputeSurfaceMesh();

  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
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
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), Material{1},
                                   Material{0});

  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
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
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), Material{1},
                                   Material{0});
  rve->ComputeSurfaceMesh();

  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
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

TEST(TestHomogenization, TestAssembleStiffnessMatrix) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();

  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);
  const MatrixXr K = homogenization->AssembleStiffnessMatrix(
      3000, unique_dof, hexahedron.at(0), hexahedron.at(1));

  ASSERT_TRUE(IsApprox(K(0, 0), 44.4445, 0.0001));
  ASSERT_TRUE(K(0, 1) < 0.0001);
}

TEST(TestHomogenization, TestAssembleLoadMatrix) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);
  const MatrixXr F = homogenization->AssembleLoadMatrix(
      1000, 3000, unique_dof, hexahedron.at(2), hexahedron.at(3));

  for (int row = 0; row < F.rows(); ++row) {
    for (int col = 0; col < F.cols(); ++col) {
      // Whole thing ends up being basically 0
      ASSERT_TRUE(F(row, col) < 0.00001);
    }
  }
}

TEST(TestHomogenization, TestAssembleLoadMatrixWithVoids) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  MatrixXr surface(10, 10);
  surface.row(0) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  surface.row(1) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(2) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(3) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(4) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(5) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(6) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(7) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(8) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(9) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  Tensor3r surface_mesh = Tensor3r::Replicate(surface, 10);

  const auto homogenization =
      std::make_shared<Homogenization>(surface_mesh, rve->PrimaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);
  const MatrixXr F = homogenization->AssembleLoadMatrix(
      1000, 3000, unique_dof, hexahedron.at(2), hexahedron.at(3));

  // Whole thing ends up being basically 0 in the sum
  ASSERT_TRUE(F.sum() < 0.00001);
}

TEST(TestHomogenization, TestComputeDisplacement) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);
  const MatrixXr F = homogenization->AssembleLoadMatrix(
      1000, 3000, unique_dof, hexahedron.at(2), hexahedron.at(3));
  const MatrixXr K = homogenization->AssembleStiffnessMatrix(
      3000, unique_dof, hexahedron.at(0), hexahedron.at(1));

  const MatrixXr X =
      homogenization->ComputeDisplacement(3000, K, F, unique_dof);

  for (int row = 0; row < X.rows(); ++row) {
    for (int col = 0; col < X.cols(); ++col) {
      // Whole thing ends up being basically 0
      ASSERT_TRUE(X(row, col) < 0.00001);
    }
  }
}

TEST(TestHomogenization, TestComputeDisplacementWithVoidNodes) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  MatrixXr surface(10, 10);
  surface.row(0) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  surface.row(1) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(2) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(3) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(4) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(5) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(6) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(7) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(8) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(9) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  Tensor3r surface_mesh = Tensor3r::Replicate(surface, 10);

  const auto homogenization =
      std::make_shared<Homogenization>(surface_mesh, rve->PrimaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);

  constexpr unsigned int n_elements = 1000;
  const MatrixX<int> edof =
      homogenization->ComputeElementDegreesOfFreedom(n_elements);
  const Tensor3i unique_nodes = homogenization->ComputeUniqueNodes(n_elements);
  const MatrixX<int> unique_dof =
      homogenization->ComputeUniqueDegreesOfFreedom(edof, unique_nodes);
  const MatrixXr F = homogenization->AssembleLoadMatrix(
      1000, 3000, unique_dof, hexahedron.at(2), hexahedron.at(3));
  const MatrixXr K = homogenization->AssembleStiffnessMatrix(
      3000, unique_dof, hexahedron.at(0), hexahedron.at(1));

  const MatrixXr X =
      homogenization->ComputeDisplacement(3000, K, F, unique_dof);

  // Just compare some of the values.
  const Real c_1 = X(4, 0);
  const Real c_2 = X(4, 1);
  const Real c_3 = X(4, 2);

  ASSERT_TRUE(IsApprox(c_1, -0.325241, 0.0001));
  ASSERT_TRUE(IsApprox(c_2, -0.0561329, 0.0001));
  ASSERT_TRUE(IsApprox(c_3, -0.0953439, 0.0001));
}

TEST(TestHomogenization, TestComputeUnitStrainParameters) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  const auto homogenization = std::make_shared<Homogenization>(
      rve->SurfaceMesh(), rve->PrimaryMaterial(), rve->SecondaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);
  const Tensor3r strain_param =
      homogenization->ComputeUnitStrainParameters(1000, hexahedron);

  const MatrixXr l0 = strain_param.Layer(0);
  VectorXr row = l0.row(0);
  Real sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(3), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(6), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(15), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));

  const MatrixXr l1 = strain_param.Layer(1);
  row = l1.row(0);
  sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(7), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(10), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(19), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(22), 1, 0.0001));

  const MatrixXr l2 = strain_param.Layer(2);
  row = l2.row(0);
  sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(14), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(17), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(20), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(23), 1, 0.0001));

  const MatrixXr l3 = strain_param.Layer(3);
  row = l3.row(0);
  sum = l3.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(6), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(9), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(21), 1, 0.0001));

  const MatrixXr l4 = strain_param.Layer(4);
  row = l4.row(0);
  sum = l4.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(13), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(16), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(19), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(22), 1, 0.0001));

  const MatrixXr l5 = strain_param.Layer(5);
  row = l5.row(0);
  sum = l5.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(12), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(15), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(21), 1, 0.0001));
}

TEST(TestHomogenization, TestComputeUnitStrainParametersWithVoids) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  MatrixXr surface(10, 10);
  surface.row(0) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  surface.row(1) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(2) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(3) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(4) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(5) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(6) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(7) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(8) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(9) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  Tensor3r surface_mesh = Tensor3r::Replicate(surface, 10);

  const auto homogenization =
      std::make_shared<Homogenization>(surface_mesh, rve->PrimaryMaterial());
  ASSERT_TRUE(homogenization.get() != nullptr);
  const auto hexahedron = homogenization->ComputeHexahedron(0.5, 0.5, 0.5);
  const Tensor3r strain_param =
      homogenization->ComputeUnitStrainParameters(1000, hexahedron);

  const MatrixXr l0 = strain_param.Layer(0);
  VectorXr row = l0.row(0);
  Real sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(3), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(6), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(15), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));

  const MatrixXr l1 = strain_param.Layer(1);
  row = l1.row(0);
  sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(7), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(10), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(19), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(22), 1, 0.0001));

  const MatrixXr l2 = strain_param.Layer(2);
  row = l2.row(0);
  sum = row.sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(14), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(17), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(20), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(23), 1, 0.0001));

  const MatrixXr l3 = strain_param.Layer(3);
  row = l3.row(0);
  sum = l3.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(6), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(9), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(21), 1, 0.0001));

  const MatrixXr l4 = strain_param.Layer(4);
  row = l4.row(0);
  sum = l4.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(13), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(16), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(19), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(22), 1, 0.0001));

  const MatrixXr l5 = strain_param.Layer(5);
  row = l5.row(0);
  sum = l5.row(0).sum();
  ASSERT_TRUE(std::fabs(sum - 4) < 0.0001);
  ASSERT_TRUE(IsApprox(row(12), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(15), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(18), 1, 0.0001));
  ASSERT_TRUE(IsApprox(row(21), 1, 0.0001));
}

TEST(TestHomogenization, TestSolverStep) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  rve->Homogenize();
  const MatrixXr constitutive_tensor = rve->ConsitutiveTensor();

  Matrix6r comp;
  comp.row(0) << 30, 10, 10, 0, 0, 0;
  comp.row(1) << 10, 30, 10, 0, 0, 0;
  comp.row(2) << 10, 10, 30, 0, 0, 0;
  comp.row(3) << 0, 0, 0, 10, 0, 0;
  comp.row(4) << 0, 0, 0, 0, 10, 0;
  comp.row(5) << 0, 0, 0, 0, 0, 10;

  ASSERT_TRUE(constitutive_tensor.isApprox(comp));
}

TEST(TestHomogenization, TestSolverStepWithVoids) {
  const auto material_1 = MaterialFromLameCoefficients(1, "one", 10, 10);
  const auto material_2 = MaterialFromLameCoefficients(2, "two", 0, 0);
  auto rve = std::make_shared<Rve>(Eigen::Vector3i(10, 10, 10), material_1,
                                   material_2);
  rve->ComputeSurfaceMesh();
  ASSERT_TRUE(rve.get() != nullptr);

  MatrixXr surface(10, 10);
  surface.row(0) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  surface.row(1) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(2) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(3) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(4) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(5) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(6) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(7) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(8) << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;
  surface.row(9) << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

  Tensor3r surface_mesh = Tensor3r::Replicate(surface, 10);

  const auto homogenization =
      std::make_shared<Homogenization>(surface_mesh, rve->PrimaryMaterial());
  homogenization->Solve();
  const MatrixXr C = homogenization->Stiffness();

  Matrix6r comp;
  comp.row(0) << 5.66474, 0.48163, 1.53659, 0, 0, 0;
  comp.row(1) << 0.48163, 5.66474, 1.53659, 0, 0, 0;
  comp.row(2) << 1.53659, 1.53659, 9.7683, 0, 0, 0;
  comp.row(3) << 0, 0, 0, 0.16672, 0, 0;
  comp.row(4) << 0, 0, 0, 0, 2.15082, 0;
  comp.row(5) << 0, 0, 0, 0, 0, 2.15082;

  ASSERT_TRUE(C.isApprox(comp));
}
