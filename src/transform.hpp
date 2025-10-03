#pragma once

#include "camera.hpp"
#include "types.hpp"
#include <concepts>

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> translation(T x, T y, T z) noexcept
{
    Mat4x4<T> t = Mat4x4<T>::identity();
    t[0, 3] = x;
    t[1, 3] = y;
    t[2, 3] = z;
    return t;
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> translation(const Vec3<T> &v) noexcept
{
    return translation(v.x(), v.y(), v.z());
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> rotation_x(T degrees) noexcept
{
    T rad = radians(degrees);
    T c = std::cos(rad);
    T s = std::sin(rad);
    Mat4x4<T> r = Mat4x4<T>::identity();
    r[1, 1] = c;
    r[1, 2] = -s;
    r[2, 1] = s;
    r[2, 2] = c;
    return r;
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> rotation_y(T degrees) noexcept
{
    T rad = radians(degrees);
    T c = std::cos(rad);
    T s = std::sin(rad);
    Mat4x4<T> r = Mat4x4<T>::identity();
    r[0, 0] = c;
    r[0, 2] = s;
    r[2, 0] = -s;
    r[2, 2] = c;
    return r;
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> rotation_z(T degrees) noexcept
{
    T rad = radians(degrees);
    T c = std::cos(rad);
    T s = std::sin(rad);
    Mat4x4<T> r = Mat4x4<T>::identity();
    r[0, 0] = c;
    r[0, 1] = -s;
    r[1, 0] = s;
    r[1, 1] = c;
    return r;
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Mat4x4<T> rotation(T xdegrees, T ydegrees, T zdegrees) noexcept
{
    const auto rx = rotation_x(xdegrees);
    const auto ry = rotation_y(ydegrees);
    const auto rz = rotation_z(zdegrees);
    return rz * ry * rx;
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> rotation(const Vec3<T> &degrees) noexcept
{
    return rotation(degrees.x(), degrees.y(), degrees.z());
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> scale(T x, T y, T z) noexcept
{
    return {
        x, T{0}, T{0}, T{0}, T{0}, y, T{0}, T{0}, T{0}, T{0}, z, T{0}, T{0}, T{0}, T{0}, T{1},
    };
}

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> scale(const Vec3<T> &v) noexcept
{
    return scale(v.x(), v.y(), v.z());
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Mat4x4<T> look_at(const Vec3<T> &eye, const Vec3<T> &target, const Vec3<T> &up) noexcept
{
    const auto forward = (target - eye).normalized().value();
    const auto right = forward.cross(up).normalized().value();
    const auto true_up = right.cross(forward);
    Mat4x4<T> m{};
    m[0, 0] = right.x();
    m[0, 1] = right.y();
    m[0, 2] = right.z();
    m[0, 3] = -right.dot(eye);

    m[1, 0] = true_up.x();
    m[1, 1] = true_up.y();
    m[1, 2] = true_up.z();
    m[1, 3] = -true_up.dot(eye);

    m[2, 0] = -forward.x();
    m[2, 1] = -forward.y();
    m[2, 2] = -forward.z();
    m[2, 3] = forward.dot(eye);

    m[3, 3] = T{1};
    return m;
}

template <std::floating_point T>
[[nodiscard]] constexpr inline Mat4x4<T> model_matrix(const Vec3<T> &position, const Vec3<T> &rotation_degrees,
                                                      const Vec3<T> &scale_) noexcept
{
    const auto t = translation(position);
    const auto r = rotation(rotation_degrees);
    const auto s = scale(scale_);
    return t * r * s;
}
