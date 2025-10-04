#pragma once

#include "check.hpp"
#include "image.hpp"
#include "types.hpp"

#include <cmath>
#include <concepts>
#include <numbers>

template <std::floating_point T> struct Camera
{
    using Degrees = T;

    Resolution resolution;
    Degrees vertical_field_of_view;
    T near_plane;
    T far_plane;
    T focal_length = T{1};
};

template <std::floating_point T> [[nodiscard]] constexpr inline T radians(T degrees) noexcept
{
    return std::numbers::pi_v<T> * degrees / T{180};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Mat4x4<T> projection_matrix(const Camera<T> &camera) noexcept
{
    const T tol{1e-6f};
    CHECK(camera.resolution.width > 0);
    CHECK(camera.far_plane > camera.near_plane + tol);

    const auto t = std::tan(radians(camera.vertical_field_of_view) / T{2});
    CHECK(std::abs(t) > tol);
    const auto sy = T{1} / t;
    const auto sx = sy * static_cast<T>(camera.resolution.height) / static_cast<T>(camera.resolution.width);

    Mat4x4<T> m{};
    m[0, 0] = sx;
    m[1, 1] = sy;
    m[2, 2] = camera.far_plane / (camera.near_plane - camera.far_plane);
    m[2, 3] = camera.far_plane * camera.near_plane / (camera.near_plane - camera.far_plane);
    m[3, 2] = -T{1};
    return m;
}
