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

auto linear_algebra::LinSpace(Real start, Real stop, unsigned int num)
    -> VectorXr {
    VectorXr interval;
    if (num == 1) {
        interval.resize(1);
        interval << (stop - start / 2) + start;
        return interval;
    }

    assert(stop > start && "STOP CANNOT BE GREATER THAN START");

    const Real div = num - 1;

    assert(div > 0 && "NUM MUST BE GREATER THAN 1");

    const Real delta = stop - start;
    interval.resize(num);

    // Initialize
    for (int i = 0; i < num; ++i) {
        interval(i) = i;
    }

    const Real step = delta / div;

    assert(step > 0 && "STEP FUNCTION INVALID CHECK DELTA AND DIV");

    interval *= step;
    interval += (VectorXr::Ones(interval.rows()) * start);

    interval(interval.rows() - 1) = stop;

    return interval;
}
