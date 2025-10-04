#pragma once
#include "types.hpp"

#include <cmath>
#include <concepts>

template <std::floating_point T> [[nodiscard]] constexpr inline size_t floor_to_size(T x) noexcept
{
    return static_cast<size_t>(std::floor(x));
}

template <std::floating_point T> [[nodiscard]] constexpr inline size_t ceil_to_size(T x) noexcept
{
    return static_cast<size_t>(std::ceil(x));
}

template <std::floating_point T>
[[nodiscard]] constexpr inline T edge_function(const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}
