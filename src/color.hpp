#pragma once

#include <array>
#include <cmath>
#include <concepts>

#include "check.hpp"
#include "types.hpp"

template <typename T> struct RGBA
{
    T r;
    T g;
    T b;
    T a;
};

template <std::floating_point T> constexpr inline RGBA<T> operator+(const RGBA<T> &c1, const RGBA<T> &c2) noexcept
{
    return {c1.r + c2.r, c1.g + c2.g, c1.b + c2.b, c1.a + c2.a};
}

template <std::floating_point T> constexpr inline RGBA<T> operator*(const RGBA<T> &c1, const RGBA<T> &c2) noexcept
{
    return {c1.r * c2.r, c1.g * c2.g, c1.b * c2.b, c1.a * c2.a};
}

template <std::floating_point T> constexpr inline RGBA<T> operator*(const RGBA<T> &c, T s) noexcept
{
    return {c.r * s, c.g * s, c.b * s, c.a * s};
}

template <std::floating_point T> constexpr inline RGBA<T> operator*(T s, const RGBA<T> &c) noexcept
{
    return c * s;
}

template <std::floating_point T> constexpr inline RGBA<T> operator/(const RGBA<T> &c, T div) noexcept
{
    T s = T{1} / div;
    return c * s;
}

template <typename T> const RGBA<T> BLACK = {T{0}, T{0}, T{0}, T{255}};

template <std::floating_point T> constexpr inline T srgb_to_linear(const T &c) noexcept
{
    if (c <= T(0.04045))
    {
        return c / T(12.92);
    }
    else
    {
        return std::pow((c + T(0.055)) / T(1.055), T(2.4));
    }
}

template <std::floating_point T> constexpr inline T linear_to_srgb(const T &c) noexcept
{
    if (c <= T(0.0031308))
    {
        return c * T(12.92);
    }
    else
    {
        return T(1.055) * std::pow(c, T(1.0 / 2.4)) - T(0.055);
    }
}

template <std::floating_point T>
[[nodiscard]] constexpr inline std::array<T, 257> generate_srgb_to_linear_table() noexcept
{
    std::array<T, 257> lut;
    for (size_t c = 0; c < 257; ++c)
    {
        lut[c] = srgb_to_linear(T(c) / T{257});
    }
    return lut;
}

template <std::floating_point T>
[[nodiscard]] constexpr inline std::array<T, 1025> generate_linear_to_srgb_table() noexcept
{
    std::array<T, 1025> lut;
    for (size_t c = 0; c < 1025; ++c)
    {
        lut[c] = linear_to_srgb(T(c) / T{1025});
    }
    return lut;
}

static const std::array<f32, 257> srgb_to_linear_table_f32 = generate_srgb_to_linear_table<f32>();
static const std::array<f64, 257> srgb_to_linear_table_f64 = generate_srgb_to_linear_table<f64>();
[[nodiscard]] constexpr inline f32 srgb_to_linearf(f32 c) noexcept
{
    c = std::clamp(c * 255.5f, 0.0f, 255.0f);
    const f32 l = std::floor(c);
    const f32 t = c - l;
    const size_t i = (size_t)l;
    const f32 interpolated = (1 - t) * srgb_to_linear_table_f32[i] + t * srgb_to_linear_table_f32[i + 1];
    return interpolated;
}

[[nodiscard]] constexpr inline f64 srgb_to_linearf(f64 c) noexcept
{
    c = std::clamp(c * 255.5, 0.0, 255.0);
    const f64 l = std::floor(c);
    const f64 t = c - l;
    const size_t i = (size_t)l;
    const f64 interpolated = (1 - t) * srgb_to_linear_table_f32[i] + t * srgb_to_linear_table_f32[i + 1];
    return interpolated;
}

template <std::floating_point T> [[nodiscard]] constexpr inline T srgb_to_linear(u8 c) noexcept
{
    (void)c;
    CHECK(false && "NOT IMPLEMENTED");
}

template <> [[nodiscard]] constexpr inline f32 srgb_to_linear(u8 c) noexcept
{
    return srgb_to_linear_table_f32[c];
}

template <> [[nodiscard]] constexpr inline f64 srgb_to_linear(u8 c) noexcept
{
    return srgb_to_linear_table_f64[c];
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<T> srgb_to_linear(const RGBA<u8> &rgb) noexcept
{
    return {
        srgb_to_linear<T>(rgb.r),
        srgb_to_linear<T>(rgb.g),
        srgb_to_linear<T>(rgb.b),
        srgb_to_linear<T>(rgb.a),
    };
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<T> srgb_to_linear(const RGBA<T> &rgb) noexcept
{
    return {
        srgb_to_linearf(rgb.r),
        srgb_to_linearf(rgb.g),
        srgb_to_linearf(rgb.b),
        srgb_to_linearf(rgb.a),
    };
}

static const std::array<f32, 1025> linear_to_srgb_table_f32 = generate_linear_to_srgb_table<f32>();
static const std::array<f64, 1025> linear_to_srgb_table_f64 = generate_linear_to_srgb_table<f64>();
constexpr inline f32 linear_to_srgbf(f32 c) noexcept
{
    c = std::clamp(c * 1024.5f, 0.0f, 1024.0f);
    const f32 l = std::floor(c);
    const f32 t = c - l;
    const size_t i = (size_t)l;
    const f32 interpolated = (1 - t) * linear_to_srgb_table_f32[i] + t * linear_to_srgb_table_f32[i + 1];
    return interpolated;
}

constexpr inline f64 linear_to_srgbf(f64 c) noexcept
{
    c = std::clamp(c * 1024.5, 0.0, 1024.0);
    const f64 l = std::floor(c);
    const f64 t = c - l;
    const size_t i = (size_t)l;
    const f64 interpolated = (1 - t) * linear_to_srgb_table_f64[i] + t * linear_to_srgb_table_f64[i + 1];
    return interpolated;
}

constexpr inline u8 linear_to_srgb(f32 c) noexcept
{
    return (u8)(linear_to_srgbf(c) * 255.0f);
}

constexpr inline u8 linear_to_srgb(f64 c) noexcept
{
    return (u8)(linear_to_srgbf(c) * 255.0);
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<u8> linear_to_srgb(const RGBA<T> &linear) noexcept
{
    return {
        linear_to_srgb(linear.r),
        linear_to_srgb(linear.g),
        linear_to_srgb(linear.b),
        linear_to_srgb(linear.a),
    };
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<T> linear_to_srgbf(const RGBA<T> &linear) noexcept
{
    return {
        linear_to_srgbf(linear.r),
        linear_to_srgbf(linear.g),
        linear_to_srgbf(linear.b),
        linear_to_srgbf(linear.a),
    };
}
