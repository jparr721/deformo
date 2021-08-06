#pragma once

#include "Numerics.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <intrin.h>
#include <iostream>
#include <numeric>
#include <ratio>
#include <sstream>
#include <type_traits>
#include <vector>

template <typename T, typename dummy = void>
struct is_printable : std::false_type {};

template <typename T>
struct is_printable<
    T, std::enable_if<std::is_same_v<
           std::remove_reference_t<decltype(std::cout << std::declval<T>())>,
           std::ostream>>> : std::true_type {};

template <typename T>
inline constexpr bool is_printable_v = is_printable<T>::value;

namespace utils {
auto OpenFile(const std::string& filename) -> std::ifstream;
auto ReadFile(const std::string& path) -> std::string;
auto Split(const std::string& input, const std::string& delimiter)
    -> std::pair<std::string, std::string>;

template <typename Derived>
void MatrixToList(std::vector<std::vector<typename Derived::Scalar>>& V,
                  const Eigen::PlainObjectBase<Derived>& M) {
    V.resize(M.rows(), std::vector<typename Derived::Scalar>(M.cols()));

    for (int row = 0; row < M.rows(); ++row) {
        for (int col = 0; col < M.cols(); ++col) {
            V[row][col] = M(row, col);
        }
    }
}

template <typename Derived>
void ListToMatrix(Eigen::PlainObjectBase<Derived>& M,
                  const std::vector<std::vector<typename Derived::Scalar>>& V) {
    M.resize(V.size(), V[0].size());

    for (int row = 0; row < V.size(); ++row) {
        for (int col = 0; col < V[0].size(); ++col) {
            M(row, col) = V[row][col];
        }
    }
}

template <typename T>
auto VectorToSTLVector(const VectorX<T>& in) -> std::vector<T> {
    std::vector<T> v;
    v.reserve(in.rows());
    for (int row = 0; row < in.rows(); ++row) {
        v.push_back(in(row));
    }

    return v;
}

template <typename T>
auto STLVectorToVector(const std::vector<T>& in) -> VectorX<T> {
    VectorX<T> v(in.size());
    for (int row = 0; row < in.size(); ++row) {
        v(row) = in.at(row);
    }

    return v;
}

template <typename Derived>
void SliceEigenVector(Eigen::PlainObjectBase<Derived>& out,
                      const Eigen::DenseBase<Derived>& in, const int start,
                      const int end) {
    assert(start < end && "YOU PROVIDED AN INVALID SLICE RANGE");
    assert(start != end && "START AND END ARE THE SAME");
    assert(start < in.rows() && "START VALUE TOO LARGE");
    assert(end <= in.rows() && "END VALUE TOO LARGE");

    out.resize((end - start) + 1, 1);
    int out_idx = -1;
    for (int i = start; i <= end; ++i) {
        out(++out_idx) = in(i);
    }
}

template <typename T>
auto SortVector(VectorX<T>& out) -> void {
    auto v = VectorToSTLVector(out);
    std::sort(v.begin(), v.end());
    
    out = STLVectorToVector(v);
}

template <typename T>
auto RemoveDuplicatesFromVector(VectorX<T>& out) -> void {
    SortVector(out);
    auto v = VectorToSTLVector(out);
    v.erase(std::unique(v.begin(), v.end()), v.end());
    out = STLVectorToVector(v);
}

/*
 * Calculates the union of two lists, removing duplicates which would originate
 * from the second list
 */
template <typename Derived>
void MatrixUnion(Eigen::PlainObjectBase<Derived>& out,
                 const Eigen::PlainObjectBase<Derived>& lhs,
                 const Eigen::PlainObjectBase<Derived>& rhs) {
    assert(lhs.cols() == rhs.cols() && "INVALID MATRIX COMBINATION");
    using T = typename Derived::Scalar;

    std::vector<std::vector<T>> v_out;
    std::vector<std::vector<T>> l;
    std::vector<std::vector<T>> r;

    MatrixToList(l, lhs);
    MatrixToList(r, rhs);

    v_out.insert(v_out.end(), r.begin(), r.end());
    v_out.insert(v_out.end(), l.begin(), l.end());
    std::sort(
        v_out.begin(), v_out.end(),
        [](const std::vector<T>& left, const std::vector<T>& right) -> bool {
            return std::accumulate(left.begin(), left.end(), 0) >
                   std::accumulate(right.begin(), right.end(), 0);
        });
    v_out.erase(std::unique(v_out.begin(), v_out.end()), v_out.end());

    ListToMatrix(out, v_out);
}

template <typename Out, typename In, typename Indices>
void SliceByIndices(Eigen::PlainObjectBase<Out>& out,
                    const Eigen::DenseBase<In>& in,
                    const Eigen::DenseBase<Indices>& rows,
                    const Eigen::DenseBase<Indices>& cols) {
    assert(rows.minCoeff() >= 0 && "ROW INDEX IS LESS THAN 0");
    assert(rows.maxCoeff() <= in.rows() &&
           "ROW INDEX IS BIGGER THAN MAX SIZE OF INPUT MATRIX");
    assert(cols.minCoeff() >= 0 && "COLUMN INDEX IS LESS THAN 0");
    assert(cols.maxCoeff() <= in.cols() &&
           "COLUMN INDEX IS BIGGER THAN MAX SIZE OF INPUT MATRIX");
    out.resize(rows.size(), cols.size());

    for (int row = 0; row < rows.size(); ++row) {
        for (int col = 0; col < cols.size(); ++col) {
            out(row, col) = in(rows(row), cols(col));
        }
    }
}

template <typename Out, typename In>
void SliceVectorByRange(Eigen::PlainObjectBase<Out>& out,
                        const Eigen::DenseBase<In>& in, int start, int end) {
    assert(in.rows() > start && in.rows() >= end &&
           "SLICE RANGE OUT OF BOUNDS");
    assert(end > start && "START HAPPENS AFTER END");
    out.resize(end - start);

    int idx = 0;
    for (int i = start; i < end; ++i, ++idx) {
        out(idx) = in(i);
    }
}

template <typename Derived>
[[nodiscard]] Real ComputeTetrahedraElementVolume(
    const Eigen::PlainObjectBase<Derived>& position_one,
    const Eigen::PlainObjectBase<Derived>& position_two,
    const Eigen::PlainObjectBase<Derived>& position_three,
    const Eigen::PlainObjectBase<Derived>& position_four) {
    const Real x1 = position_one.x();
    const Real y1 = position_one.y();
    const Real z1 = position_one.z();

    const Real x2 = position_two.x();
    const Real y2 = position_two.y();
    const Real z2 = position_two.z();

    const Real x3 = position_three.x();
    const Real y3 = position_three.y();
    const Real z3 = position_three.z();

    const Real x4 = position_four.x();
    const Real y4 = position_four.y();
    const Real z4 = position_four.z();

    Eigen::Matrix4f V;
    V.row(0) << 1, x1, y1, z1;
    V.row(1) << 1, x2, y2, z2;
    V.row(2) << 1, x3, y3, z3;
    V.row(3) << 1, x4, y4, z4;

    return V.determinant() / 6;
}

template <typename Derived>
Vector2<unsigned> Shape(const Eigen::PlainObjectBase<Derived>& in) {
    return {in.rows(), in.cols()};
}

template <typename T> void GTestDebugPrint(T value, bool sep = true) {
    std::cerr << value << std::endl;
    if (sep) {
        std::cerr << "============" << std::endl;
    }
}

void FindMaxVertices(std::vector<unsigned int>& indices,
                     const VectorXr& positions);

namespace stopwatch {

/*Some performance analysis tools*/
// An implementation of the 'TrivialClock' concept using the rdtscp instruction.
struct rdtscp_clock {
    using rep = std::uint64_t;
    using period = std::ratio<1>;
    using duration = std::chrono::duration<rep, period>;
    using time_point = std::chrono::time_point<rdtscp_clock, duration>;

    static auto now() noexcept -> time_point {
#ifdef __linux__
        std::uint32_t hi, lo;
        __asm__ __volatile__("rdtscp" : "=d"(hi), "=a"(lo));
        return time_point(
            duration((static_cast<std::uint64_t>(hi) << 32) | lo));

#else
        std::uint32_t aux;
        const std::uint64_t lo = __rdtscp(&aux);
        return time_point(duration(lo));
#endif
    }
};

// A timer using the specified clock.
template <class Clock = std::chrono::system_clock> struct timer {
    using TimePoint = typename Clock::time_point;
    using Duration = typename Clock::duration;

    explicit timer(const Duration duration) noexcept
        : expiry(Clock::now() + duration) {}
    explicit timer(const TimePoint expiry) noexcept : expiry(expiry) {}

    bool done(TimePoint now = Clock::now()) const noexcept {
        return now >= expiry;
    }

    auto remaining(TimePoint now = Clock::now()) const noexcept -> Duration {
        return expiry - now;
    }

    const TimePoint expiry;
};

template <class Clock = std::chrono::system_clock>
constexpr auto make_timer(typename Clock::duration duration) -> timer<Clock> {
    return timer<Clock>(duration);
}

// Times how long it takes a function to execute using the specified clock.
template <class Clock = rdtscp_clock, class Func>
auto time(Func&& function) -> typename Clock::duration {
    const auto start = Clock::now();
    function();
    return Clock::now() - start;
}

// Samples the given function N times using the specified clock.
template <std::size_t N, class Clock = rdtscp_clock, class Func>
auto sample(Func&& function) -> std::array<typename Clock::duration, N> {
    std::array<typename Clock::duration, N> samples;

    for (std::size_t i = 0u; i < N; ++i) {
        samples[i] = time<Clock>(function);
    }

    std::sort(samples.begin(), samples.end());
    return samples;
}
} // namespace stopwatch
} // namespace utils
