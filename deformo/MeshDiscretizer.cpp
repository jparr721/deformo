#include "MeshDiscretizer.h"
#include "Utils.h"
#include <cmath>

auto MeshDiscretizer::Discretize(const Real resolution) -> void {
    assert(resolution > 0 && resolution <= 1 && "RESOLUTION TOO LARGE");

    resolution_ = resolution;

    std::vector<VectorXr> expanded_forms;

    for (int row = 0; row < volumes_.rows(); row += 3) {
        const auto r = Eigen::Vector3i(volumes_(row), volumes_(row + 1),
                                       volumes_(row + 2));
        const Vector3r pt_1 = RowSlice(mesh_, r.x());
        const Vector3r pt_2 = RowSlice(mesh_, r.y());
        const Vector3r pt_3 = RowSlice(mesh_, r.z());

        const VectorXr x_y = RangeLerp(pt_1, pt_2);
        const VectorXr y_z = RangeLerp(pt_2, pt_3);
        const VectorXr x_z = RangeLerp(pt_1, pt_3);

        expanded_forms.emplace_back(x_y);
        expanded_forms.emplace_back(y_z);
        expanded_forms.emplace_back(x_z);
    }

    discretized_mesh_.resize(3 * expanded_forms.size());

    const long long expanded_forms_size = expanded_forms[0].rows();
    assert(expanded_forms_size > 0);
    for (long long segment = 0l; segment < expanded_forms.size(); ++segment) {
        const VectorXr expanded_form = expanded_forms.at(segment);
        discretized_mesh_.segment(segment * expanded_forms_size,
                                  expanded_forms_size)
            << expanded_form;

        // TODO(@jparr721) Volumes.
        discretized_volumes_ = volumes_; // Replace this
    }
}

auto MeshDiscretizer::ToMesh() -> std::shared_ptr<Mesh> {
    return std::make_shared<Mesh>(discretized_mesh_, discretized_volumes_);
}

auto MeshDiscretizer::RangeLerp(const Vector3r& a, const Vector3r& b) const
    -> VectorXr {
    const auto lerp_at_resolution = [](const Vector3r& aa, const Vector3r& bb,
                                       const Real resolution) -> VectorXr {
        const int output_size = std::log10(resolution) * 3.f;

        int segment = 0;
        Real i = 0;
        VectorXr output;
        output.resize(output_size);

        while (i < static_cast<Real>(1.1)) {
            output.segment(segment, 3) << linear_algebra::Lerp(aa, bb, i);
            i += resolution;
            segment += 3;
        }

        return output;
    };

    if (a.sum() > b.sum()) {
        return lerp_at_resolution(b, a, resolution_);
    }

    return lerp_at_resolution(a, b, resolution_);
}

auto MeshDiscretizer::RowSlice(const VectorXr& data, int index) -> VectorXr {
    index = index * 3;
    VectorXr point;
    utils::SliceVectorByRange(point, data, index, index + 2);
    return point;
}
