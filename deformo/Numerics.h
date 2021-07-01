#pragma once

#include <vector>

namespace utils::numerics {
auto OneDimensionalLinearInterpolation(float low, float high, float interval)
    -> std::vector<float>;
} // namespace utils::numerics
