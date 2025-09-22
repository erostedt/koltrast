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

template <typename T> using Mat4x4 = Matrix<T, 4, 4>;
using Mat4x4f = Mat4x4<f32>;
