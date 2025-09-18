#pragma once
#include "matrix.hpp"
#include <concepts>

template <typename T> using Mat4x4 = Matrix<T, 4, 4>;
template <typename T> using Vec3 = Vector<T, 3>;

template <std::floating_point T> [[nodiscard]] constexpr inline Mat4x4<T> translation(T x, T y, T z) noexcept
{
    Mat4x4<T> t = Mat4x4<T>::identity();
    t[0, 3] = x;
    t[1, 3] = y;
    t[2, 3] = z;
    return t;
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
