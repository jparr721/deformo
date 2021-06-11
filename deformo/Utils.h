#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <intrin.h>
#include <iostream>
#include <numeric>
#include <ratio>
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
/*
 * Calculates the union of two lists, removing duplicates which would originate
 * from the second list
 */
template <typename Derived>
void MatrixUnion(Eigen::PlainObjectBase<Derived>& out,
                 const Eigen::PlainObjectBase<Derived>& lhs,
                 const Eigen::PlainObjectBase<Derived>& rhs) {
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

template <typename T> void GTestDebugPrint(T value) {
    std::cerr << value << std::endl;
}

void FindMaxVertices(std::vector<unsigned int>& indices,
                     const Eigen::VectorXf& positions);

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
