#pragma once

#include "matrix.hpp"

#include <cstdint>

using u8 = uint8_t;
using u32 = uint32_t;
using u64 = uint64_t;

using i32 = int32_t;

using f32 = float;
using f64 = double;

template <typename T> using Vec2 = Vector<T, 2>;
template <typename T> using Vec3 = Vector<T, 3>;
template <typename T> using Vec4 = Vector<T, 4>;

using Vec2f = Vec2<f32>;
using Vec3f = Vec3<f32>;
using Vec4f = Vec4<f32>;

using Vec2d = Vec2<f64>;
using Vec3d = Vec3<f64>;
using Vec4d = Vec4<f64>;

template <typename T> using Mat4x4 = Matrix<T, 4, 4>;
template <typename T> using Mat3x3 = Matrix<T, 3, 3>;
using Mat4x4f = Mat4x4<f32>;
using Mat3x3f = Mat3x3<f32>;

using Mat4x4d = Mat4x4<f64>;
using Mat3x3d = Mat3x3<f64>;
