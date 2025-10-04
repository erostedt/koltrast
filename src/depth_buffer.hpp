#pragma once

#include <concepts>
#include <execution>

#include "image.hpp"

template <std::floating_point T, size_t AARows = 1, size_t AACols = AARows>
    requires(AARows > 0) && (AACols > 0)
using DepthBuffer = Image<Matrix<T, AARows, AACols>>;

template <std::floating_point T, size_t AARows = 1, size_t AACols = AARows>
constexpr inline void reset_depth_buffer(DepthBuffer<T, AARows, AACols> &buffer) noexcept
{
    using namespace std;
    Matrix<T, AARows, AACols> invalid;
    fill(begin(invalid), end(invalid), numeric_limits<T>::infinity());
    fill(execution::par_unseq, begin(buffer), end(buffer), invalid);
}

template <std::floating_point T, size_t AARows = 1, size_t AACols = AARows>
[[nodiscard]] inline DepthBuffer<T, AARows, AACols> create_depth_buffer(size_t width, size_t height)
{
    using namespace std;
    DepthBuffer<T, AARows, AACols> buffer(width, height);
    reset_depth_buffer(buffer);
    return buffer;
}
