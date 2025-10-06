#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <functional>

#include "check.hpp"
#include "matrix.hpp"
#include "types.hpp"

template <typename T, size_t Channels>
    requires(Channels > 0)
struct Color : public Array<T, Channels>
{
    using Array<T, Channels>::elements;
    using Array<T, Channels>::Array;
    using Array<T, Channels>::operator[];

    [[nodiscard]] constexpr inline T &r() noexcept
    {
        return elements[0];
    }

    [[nodiscard]] constexpr inline const T &r() const noexcept
    {
        return elements[0];
    }

    [[nodiscard]] constexpr inline T &g() noexcept
        requires(Channels > 1)
    {
        return elements[1];
    }

    [[nodiscard]] constexpr inline const T &g() const noexcept
        requires(Channels > 1)
    {
        return elements[1];
    }

    [[nodiscard]] constexpr inline T &b() noexcept
        requires(Channels > 2)
    {
        return elements[2];
    }

    [[nodiscard]] constexpr inline const T &b() const noexcept
        requires(Channels > 2)
    {
        return elements[2];
    }

    [[nodiscard]] constexpr inline T &a() noexcept
        requires(Channels > 3)
    {
        return elements[3];
    }

    [[nodiscard]] constexpr inline const T &a() const noexcept
        requires(Channels > 3)
    {
        return elements[3];
    }
};

template <std::floating_point T, size_t Channels>
[[nodiscard]] constexpr inline Color<T, Channels> operator*(const Color<T, Channels> &left,
                                                            const Color<T, Channels> &right) noexcept
{
    return map(left, right, std::multiplies<T>{});
}

template <typename T> using RGB = Color<T, 3>;
template <typename T> using RGBA = Color<T, 4>;

template <typename T> const RGBA<T> BLACK = {T{0}, T{0}, T{0}, T{1}};

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
        srgb_to_linear<T>(rgb.r()),
        srgb_to_linear<T>(rgb.g()),
        srgb_to_linear<T>(rgb.b()),
        srgb_to_linear<T>(rgb.a()),
    };
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<T> srgb_to_linear(const RGBA<T> &rgb) noexcept
{
    return {
        srgb_to_linearf(rgb.r()),
        srgb_to_linearf(rgb.g()),
        srgb_to_linearf(rgb.b()),
        srgb_to_linearf(rgb.a()),
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
        linear_to_srgb(linear.r()),
        linear_to_srgb(linear.g()),
        linear_to_srgb(linear.b()),
        linear_to_srgb(linear.a()),
    };
}

template <std::floating_point T> [[nodiscard]] constexpr inline RGBA<T> linear_to_srgbf(const RGBA<T> &linear) noexcept
{
    return {
        linear_to_srgbf(linear.r()),
        linear_to_srgbf(linear.g()),
        linear_to_srgbf(linear.b()),
        linear_to_srgbf(linear.a()),
    };
}
