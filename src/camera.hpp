#pragma once

#include "check.hpp"
#include "image.hpp"
#include "matrix.hpp"
#include "types.hpp"

#include <cmath>
#include <concepts>
#include <execution>
#include <numbers>
#include <vector>

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

template <std::floating_point T>
constexpr inline void model_to_world(const std::vector<Vec3<T>> &model_vertices, const Mat4x4<T> &model_matrix,
                                     std::vector<Vec4<T>> &world_vertices) noexcept
{
    world_vertices.resize(model_vertices.size());

    using namespace std;
    transform(execution::par_unseq, begin(model_vertices), end(model_vertices), begin(world_vertices),
              [&](const Vec3<T> &vertex) { return model_matrix * Vec4<T>{vertex.x(), vertex.y(), vertex.z(), T{1}}; });
}

template <std::floating_point T>
constexpr inline void model_to_world(const std::vector<Vec3<T>> &model_normals, const Mat4x4<T> &model_matrix,
                                     std::vector<Vec3<T>> &world_normals) noexcept
{
    Mat3x3<T> mat3 = {model_matrix[0, 0], model_matrix[0, 1], model_matrix[0, 2],
                      model_matrix[1, 0], model_matrix[1, 1], model_matrix[1, 2],
                      model_matrix[2, 0], model_matrix[2, 1], model_matrix[2, 2]};

    Mat3x3<T> normal_matrix = mat3.inverse().value().transposed();

    world_normals.resize(model_normals.size());

    using namespace std;
    transform(execution::par_unseq, begin(model_normals), end(model_normals), begin(world_normals),
              [&](const Vec3<T> &normal) { return *(normal_matrix * normal).normalized(); });
}

template <std::floating_point T>
void world_to_screen(const std::vector<Vec4<T>> &world_vertices, const Mat4x4<T> &vp, const Resolution &resolution,
                     std::vector<Vec4<T>> &screen_vertices)
{
    using namespace std;

    auto project_single = [&](const Vec4<T> &world_point) -> Vec4<T> {
        const auto clip = vp * world_point;
        CHECK(abs(clip.w()) > T(1e-6));
        const Vec3<T> ndc = {clip.x() / clip.w(), clip.y() / clip.w(), clip.z() / clip.w()};
        const T sx = (ndc.x() * T{0.5} + T{0.5}) * static_cast<T>(resolution.width);
        const T sy = (T{1} - (ndc.y() * T{0.5} + T{0.5})) * static_cast<T>(resolution.height);
        return {sx, sy, ndc.z(), clip.w()};
    };

    screen_vertices.resize(world_vertices.size());
    transform(execution::par_unseq, begin(world_vertices), end(world_vertices), begin(screen_vertices), project_single);
}
