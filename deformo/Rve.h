#pragma once

#include "DeformoAssert.h"
#include "Homogenization.h"
#include "ImplicitSurfaceGenerator.h"
#include "MarchingCubes.h"
#include "Material.h"
#include "Numerics.h"
#include <memory>
#include <string>
#include <utility>

class Rve {
  public:
    Rve(const Vector3<int>& size, const Material& material_1);
    Rve(const Vector3<int>& size, const Material& material_1,
        const Material& material_2);

    auto Homogenize() -> void;
    auto ComputeRenderableMesh(MatrixXr& V, MatrixX<int>& F) -> void;
    auto ComputeSurfaceMesh() -> void;
    auto ComputeSurfaceMesh(const Vector3<int>& inclusion_size,
                            const int n_inclusions, const bool is_isotropic)
        -> void;

    auto Height() const noexcept -> int { return height_; }
    auto Width() const noexcept -> int { return width_; }
    auto Depth() const noexcept -> int { return depth_; }
    auto ConsitutiveTensor() const -> Matrix6r { return C_; }
    auto SurfaceMesh() const -> Tensor3r { return surface_mesh_; }
    auto PrimaryMaterial() const -> Material { return material_1_; }
    auto SecondaryMaterial() const -> Material { return material_2_; }
    auto Homogenized() const noexcept -> const std::unique_ptr<Homogenization>& {
        return homogenization_;
    }

  private:
    static constexpr unsigned int kCellLength = 1;

    bool is_homogenized_ = false;
    bool contains_surface_mesh_ = false;

    unsigned int height_ = 0;
    unsigned int width_ = 0;
    unsigned int depth_ = 0;

    DeformoAssertion deformo_assertion_;

    Material material_1_;
    Material material_2_;

    std::unique_ptr<ImplicitSurfaceGenerator<Real>> generator_;
    std::unique_ptr<MarchingCubes> marching_cubes_;
    std::unique_ptr<Homogenization> homogenization_;

    Matrix6r C_;
    Tensor3r surface_mesh_;
};
