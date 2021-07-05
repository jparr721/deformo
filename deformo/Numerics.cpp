#include "Numerics.h"

#include <cassert>
#include <vector>

auto linear_algebra::OneDimensionalLinearInterpolation(const Real low,
                                                       const Real high,
                                                       const Real interval)
    -> std::vector<Real> {
    assert(high - low >= 1 &&
           "Interpolation Not Supported For Values of Diff < 1");
    std::vector<Real> values;

    for (Real i = low; i <= high;) {
        values.push_back(i);
        i += interval;
    }

    // If the interval is wrong we may overshoot the highest value, add it if
    // so.
    if (values.at(values.size() - 1) < high) {
        values.push_back(high);
    }

    return values;
}
