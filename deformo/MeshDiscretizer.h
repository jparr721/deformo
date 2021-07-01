#pragma once

#include "Numerics.h"
#include <Eigen/Dense>
#include <vector>

template <typename DerivedMesh, typename DerivedVolumes> class MeshDiscretizer {
  public:
    MeshDiscretizer(
        const Eigen::PlainObjectBase<DerivedMesh> tethrahedral_mesh,
        const Eigen::PlainObjectBase<DerivedVolumes> tetrahedral_volumes)
        : mesh_(tethrahedral_mesh), volumes_(tetrahedral_volumes) {}

    constexpr auto Reoslution() const noexcept -> Real { return resolution_; }

    auto Discretize(const Real resolution)
        -> Eigen::PlainObjectBase<DerivedMesh> {
        assert(resolution > 0 && resolution <= 1 && "RESOLUTION TOO LARGE");

        resolution_ = resolution;

        return {};
    }

  private:
    Real resolution_;

    Eigen::PlainObjectBase<DerivedMesh> mesh_;
    Eigen::PlainObjectBase<DerivedVolumes> volumes_;
    Eigen::PlainObjectBase<DerivedMesh> discretized_mesh_;

    auto RangeLerp(const Real resolution) -> std::vector<Real> {}
};
