#pragma once

#include "check.hpp"
#include "matrix.hpp"
#include "types.hpp"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <execution>
#include <numbers>
#include <vector>

template <typename T> struct RGB
{
    T r;
    T g;
    T b;
};

template <typename T> const RGB<T> BLACK = {T{0}, T{0}, T{0}};

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
    T focal_length = T{1};
};

template <typename T> class Image
{
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

  public:
    Image(size_t width, size_t height) : res(width, height), pixels(width * height)
    {
    }

    [[nodiscard]] constexpr inline size_t width() const noexcept
    {
        return res.width;
    }

    [[nodiscard]] constexpr inline size_t height() const noexcept
    {
        return res.height;
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return width() * height();
    }

    [[nodiscard]] constexpr inline const Resolution &resolution() const noexcept
    {
        return resolution;
    }

    [[nodiscard]] constexpr inline T &operator[](size_t i) noexcept
    {
        return pixels[i];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t x, size_t y) noexcept
    {
        return operator[](y * width() + x);
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t i) const
    {
        return pixels[i];
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t x, size_t y) const noexcept
    {
        return operator[](y * width() + x);
    }

    [[nodiscard]] constexpr inline iterator begin() noexcept
    {
        return std::begin(pixels);
    }

    [[nodiscard]] constexpr inline iterator end() noexcept
    {
        return std::end(pixels);
    }

    [[nodiscard]] constexpr inline const_iterator begin() const noexcept
    {
        return std::cbegin(pixels);
    }

    [[nodiscard]] constexpr inline const_iterator end() const noexcept
    {
        return std::cend(pixels);
    }

  private:
    Resolution res;
    std::vector<T> pixels;
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
constexpr inline void model_to_world(const std::vector<Vec3<T>> &model_vertices,
                                     const std::vector<Vec3<T>> &model_normals, const Mat4x4<T> &model_matrix,
                                     std::vector<Vec4<T>> &world_vertices, std::vector<Vec3<T>> &world_normals) noexcept
{
    Mat3x3<T> mat3 = {model_matrix[0, 0], model_matrix[0, 1], model_matrix[0, 2],
                      model_matrix[1, 0], model_matrix[1, 1], model_matrix[1, 2],
                      model_matrix[2, 0], model_matrix[2, 1], model_matrix[2, 2]};

    Mat3x3<T> normal_matrix = mat3.inverse().value().transposed();

    world_vertices.resize(model_vertices.size());
    world_normals.resize(model_normals.size());

    using namespace std;
    transform(execution::par_unseq, begin(model_vertices), end(model_vertices), begin(world_vertices),
              [&](const Vec3<T> &vertex) { return model_matrix * Vec4<T>{vertex.x(), vertex.y(), vertex.z(), T{1}}; });

    transform(execution::par_unseq, begin(model_normals), end(model_normals), begin(world_normals),
              [&](const Vec3<T> &normal) { return *(normal_matrix * normal).normalized(); });
}

template <std::floating_point T>
Vec4<T> project_to_screen(const Vec4<T> &world_point, const Mat4x4<T> &vp, const Resolution &res)
{
    const auto clip = vp * world_point;
    CHECK(std::abs(clip.w()) > T(1e-6));
    const Vec3<T> ndc = {clip.x() / clip.w(), clip.y() / clip.w(), clip.z() / clip.w()};
    const T sx = (ndc.x() * T{0.5} + T{0.5}) * static_cast<T>(res.width);
    const T sy = (T{1} - (ndc.y() * T{0.5} + T{0.5})) * static_cast<T>(res.height);

    return {sx, sy, ndc.z(), clip.w()};
}

template <std::floating_point T>
void project_to_screen(const std::vector<Vec4<T>> &world_vertices, const Mat4x4<T> &vp, const Resolution &resolution,
                       std::vector<Vec4<T>> &screen_vertices)
{
    using namespace std;
    screen_vertices.resize(world_vertices.size());
    transform(execution::par_unseq, begin(world_vertices), end(world_vertices), begin(screen_vertices),
              [&](const Vec4<T> &vertex) { return project_to_screen(vertex, vp, resolution); });
}
