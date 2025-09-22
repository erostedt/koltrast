#pragma once
#include "types.hpp"

#include <cmath>

[[nodiscard]] constexpr inline size_t floor_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::floor(x));
}

[[nodiscard]] constexpr inline size_t ceil_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::ceil(x));
}

[[nodiscard]] constexpr inline f32 edge_function(Vec4f p1, Vec4f p2, Vec4f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}
