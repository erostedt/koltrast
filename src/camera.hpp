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

template <std::floating_point T> struct ViewPort
{
    Vec3<T> pixel_delta_u;
    Vec3<T> pixel_delta_v;
    Vec3<T> pixel_top_left;
};

template <std::floating_point T>
ViewPort<T> view_port(const Camera<T> &camera, const Vec3<T> &position, const Mat4x4<T> &view)
{
    T aspect_ratio = (T)camera.resolution.width / (T)camera.resolution.height;

    T angle = radians(camera.vertical_field_of_view);
    T height = std::tan(angle / T{2});
    T viewport_height = T{2} * height * camera.focal_length;
    T viewport_width = viewport_height * aspect_ratio;

    Vec3<T> right = {view[0, 0], view[0, 1], view[0, 2]};
    Vec3<T> up = {view[1, 0], view[1, 1], view[1, 2]};
    Vec3<T> forward = {view[2, 0], view[2, 1], view[2, 2]};

    Vec3<T> viewport_u = viewport_width * right;
    Vec3<T> viewport_v = -viewport_height * up;

    Vec3<T> viewport_center = position - camera.focal_length * forward;
    Vec3<T> viewport_upper_left = viewport_center - viewport_u / T{2} - viewport_v / T{2};

    Vec3<T> pixel_delta_u = viewport_u / (T)camera.resolution.width;
    Vec3<T> pixel_delta_v = viewport_v / (T)camera.resolution.height;

    Vec3<T> pixel_top_left = viewport_upper_left + T(0.5) * (pixel_delta_u + pixel_delta_v);
    return {pixel_delta_u, pixel_delta_v, pixel_top_left};
}
