#pragma once

#include "matrix.hpp"
#include "types.hpp"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

struct RGB
{
    u8 r;
    u8 g;
    u8 b;
};

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

template <typename T> class Image
{
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

  public:
    Image(size_t width, size_t height) : resolution(width, height), pixels(width * height)
    {
    }

    [[nodiscard]] constexpr inline size_t width() const noexcept
    {
        return resolution.width;
    }

    [[nodiscard]] constexpr inline size_t height() const noexcept
    {
        return resolution.height;
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
    Resolution resolution;
    std::vector<T> pixels;
};

template <std::floating_point T> [[nodiscard]] constexpr inline T radians(T degrees) noexcept
{
    return std::numbers::pi_v<T> * degrees / T{180};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Matrix<T, 4, 4> projection_matrix(const Camera<T> &camera) noexcept
{
    const T tol{1e-6f};
    CHECK(camera.resolution.width > 0);
    CHECK(camera.far_plane > camera.near_plane + tol);

    const auto t = std::tan(radians(camera.vertical_field_of_view) / T{2});
    CHECK(std::abs(t) > tol);
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
    CHECK(std::abs(clip.w()) > T{1e-6f});
    const Vector<T, 3> ndc = {clip.x() / clip.w(), clip.y() / clip.w(), clip.z() / clip.w()};
    const T sx = (ndc.x() * T{0.5} + T{0.5}) * static_cast<T>(res.width);
    const T sy = (T{1} - (ndc.y() * T{0.5} + T{0.5})) * static_cast<T>(res.height);

    return {sx, sy, ndc.z()};
}
