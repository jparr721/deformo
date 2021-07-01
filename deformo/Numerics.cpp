#include "Numerics.h"

#include <cassert>
#include <vector>

auto utils::numerics::OneDimensionalLinearInterpolation(const float low,
                                                        const float high,
                                                        const float interval)
    -> std::vector<float> {
    assert(high - low >= 1 &&
           "Interpolation Not Supported For Values of Diff < 1");
    std::vector<float> values;

    for (float i = low; i <= high;) {
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
