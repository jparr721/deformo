#pragma once

#include "Mesh.h"
#include "Numerics.h"
#include <utility>

class MeshDiscretizer {
  public:
    MeshDiscretizer(VectorXr tethrahedral_mesh,
                    Eigen::VectorXi tetrahedral_volumes)
        : mesh_(std::move(tethrahedral_mesh)),
          volumes_(std::move(tetrahedral_volumes)) {}

    auto Discretize(Real resolution = static_cast<Real>(0.1)) -> void;
    auto ToMesh() -> std::shared_ptr<Mesh>;
    auto ToExpandedForm(const VectorXr& data_dim) -> Eigen::Tensor<Real, 4>;

    [[nodiscard]] constexpr auto Resolution() const noexcept -> Real {
        return resolution_;
    }

  private:
    Real resolution_;

    VectorXr mesh_;
    Eigen::VectorXi volumes_;

    VectorXr discretized_mesh_;
    Eigen::VectorXi discretized_volumes_;

    auto RangeLerp(const Vector3r& a, const Vector3r& b) const -> VectorXr;

    auto RowSlice(const VectorXr& data, int index) -> VectorXr;
};
