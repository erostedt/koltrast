#pragma once

#include "matrix.hpp"
#include <concepts>
#include <cstddef>

struct Resolution
{
    size_t width = 640;
    size_t height = 480;
};

template <std::floating_point T> struct Camera
{
    using Degrees = T;

    Resolution resolution;
    Degrees vertical_field_of_view;
    T near_plane;
    T far_plane;
};

template <std::floating_point T> [[nodiscard]] constexpr inline T radians(T degrees) noexcept
{
    return std::numbers::pi_v<T> * degrees / T{180};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Matrix<T, 4, 4> projection_matrix(const Camera<T> &camera) noexcept
{
    const T tol{1e-6f};
    assert(camera.resolution.width > 0 && "No width");
    assert(camera.far_plane > camera.near_plane + tol && "Near and far plane are the same");

    const auto t = std::tan(radians(camera.vertical_field_of_view) / T{2});
    assert(std::abs(t) > tol && "Bad vfov");
    const auto sy = T{1} / t;
    const auto sx = sy * static_cast<T>(camera.resolution.height) / static_cast<T>(camera.resolution.width);

    Matrix<T, 4, 4> m{};
    m[0, 0] = sx;
    m[1, 1] = sy;
    m[2, 2] = camera.far_plane / (camera.near_plane - camera.far_plane);
    m[2, 3] = camera.far_plane * camera.near_plane / (camera.near_plane - camera.far_plane);
    m[3, 2] = -T{1};
    return m;
}

template <std::floating_point T>
Vector<T, 3> project_to_screen(const Vector<T, 3> &world_point, const Matrix<T, 4, 4> &view,
                               const Matrix<T, 4, 4> &proj, const Resolution &res)
{
    const Vector<T, 4> world_point_homo = {world_point.x(), world_point.y(), world_point.z(), T{1}};
    const auto camera_point_homo = view * world_point_homo;
    const Vector<T, 4> clip = proj * camera_point_homo;
    assert(std::abs(clip.w()) > T{1e-6f});
    const Vector<T, 3> ndc = {clip.x() / clip.w(), clip.y() / clip.w(), clip.z() / clip.w()};
    const T sx = (ndc.x() * T{0.5} + T{0.5}) * static_cast<T>(res.width);
    const T sy = (T{1} - (ndc.y() * T{0.5} + T{0.5})) * static_cast<T>(res.height);

    return {sx, sy, ndc.z()};
}
