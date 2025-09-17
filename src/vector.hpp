#pragma once

#include "elementwise.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <numeric>

template <typename T, size_t Rows>
    requires std::is_signed_v<T>
class Vector
{
    using value_type = typename std::array<T, Rows>::value_type;
    using pointer = typename std::array<T, Rows>::pointer;
    using const_pointer = typename std::array<T, Rows>::const_pointer;
    using reference = typename std::array<T, Rows>::reference;
    using const_reference = typename std::array<T, Rows>::const_reference;
    using iterator = typename std::array<T, Rows>::iterator;
    using const_iterator = typename std::array<T, Rows>::const_iterator;
    using size_type = typename std::array<T, Rows>::size_type;
    using difference_type = typename std::array<T, Rows>::difference_type;
    using reverse_iterator = typename std::array<T, Rows>::reverse_iterator;
    using const_reverse_iterator = typename std::array<T, Rows>::const_reverse_iterator;

  public:
    constexpr Vector() : elements{}
    {
    }

    template <typename... U> constexpr Vector(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == Rows, "Wrong number of elements");
        static_assert((std::is_same_v<U, T> && ...), "All arguments must be exactly T");
    }

    constexpr Vector(const std::array<T, Rows> &arr) : elements(arr)
    {
    }

    [[nodiscard]] constexpr inline iterator begin() noexcept
    {
        return std::begin(elements);
    }

    [[nodiscard]] constexpr inline const_iterator begin() const noexcept
    {
        return std::begin(elements);
    }

    [[nodiscard]] constexpr inline iterator end() noexcept
    {
        return std::end(elements);
    }

    [[nodiscard]] constexpr inline const_iterator end() const noexcept
    {
        return std::end(elements);
    }

    [[nodiscard]] constexpr inline bool operator==(const Vector &v) const noexcept
        requires std::equality_comparable<T>
    {
        return elements == v.elements;
    }

    [[nodiscard]] constexpr inline bool operator!=(const Vector &v) const noexcept
        requires std::equality_comparable<T>
    {
        return !(*this == v);
    }

    [[nodiscard]] static constexpr inline size_t size() noexcept
    {
        return Rows;
    }

    [[nodiscard]] constexpr inline T dot(const Vector &v) const noexcept
        requires Multiplicable<T, T> && Addable<T, T>
    {
        return std::inner_product(this->begin(), this->end(), std::begin(v), T{});
    }

    [[nodiscard]] constexpr inline T squared_length() const noexcept
        requires std::floating_point<T>
    {
        return dot(*this);
    }

    [[nodiscard]] constexpr inline T length() const noexcept
        requires std::floating_point<T>
    {
        return std::sqrt(squared_length());
    }

    [[nodiscard]] constexpr inline T x() const noexcept
        requires(Rows >= 1 && Rows <= 4)
    {
        return elements[0];
    }

    [[nodiscard]] constexpr inline T y() const noexcept
        requires(Rows >= 2 && Rows <= 4)
    {
        return elements[1];
    }

    [[nodiscard]] constexpr inline T z() const noexcept
        requires(Rows >= 3 && Rows <= 4)
    {
        return elements[2];
    }

    [[nodiscard]] constexpr inline T w() const noexcept
        requires(Rows == 4)
    {
        return elements[3];
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t index) const noexcept
    {
        return elements[index];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t index) noexcept
    {
        return elements[index];
    }

    [[nodiscard]] constexpr inline Vector cross(const Vector &v) const noexcept
        requires(Rows == 3) && Multiplicable<T, T> && Subtractable<T, T>
    {
        const auto &u = *this;
        const auto x = u.y() * v.z() - u.z() * v.y();
        const auto y = u.z() * v.x() - u.x() * v.z();
        const auto z = u.x() * v.y() - u.y() * v.x();
        return {x, y, z};
    }

    [[nodiscard]] constexpr inline T squared_distance(const Vector &to) const noexcept
        requires std::floating_point<T>
    {
        return (*this - to).squared_length();
    }

    [[nodiscard]] constexpr inline T distance(const Vector &to) const noexcept
        requires std::floating_point<T>
    {
        return std::sqrt(squared_distance(to));
    }

    std::array<T, Rows> elements{};
};

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator+(const Vector<T, Rows> &left,
                                                         const Vector<T, Rows> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right.elements, std::plus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator+(const Vector<T, Rows> &left, const T &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right, std::plus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator+(const T &left, const Vector<T, Rows> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left, right.elements, std::plus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator-(const Vector<T, Rows> &left,
                                                         const Vector<T, Rows> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right.elements, std::minus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator-(const Vector<T, Rows> &left, const T &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right, std::minus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator-(const T &left, const Vector<T, Rows> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left, right.elements, std::minus<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator-(const Vector<T, Rows> &v) noexcept
    requires Negatable<T>
{
    using namespace std;
    auto out = v;
    transform(begin(v), end(v), begin(out), std::negate<T>{});
    return out;
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator*(const Vector<T, Rows> &left, T right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left.elements, right, std::multiplies<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator*(T left, const Vector<T, Rows> &right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left, right.elements, std::multiplies<T>{});
}

template <typename T, size_t Rows>
[[nodiscard]] constexpr inline Vector<T, Rows> operator/(const Vector<T, Rows> &numerator, T denominator) noexcept
    requires Divisable<T, T>
{
    return Elementwise(numerator.elements, denominator, std::divides<T>{});
}
