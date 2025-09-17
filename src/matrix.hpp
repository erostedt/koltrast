#pragma once

#include "elementwise.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <functional>

template <typename T, size_t Rows, size_t Cols>
    requires std::is_signed_v<T>
class Matrix
{
    using value_type = typename std::array<T, Rows * Cols>::value_type;
    using pointer = typename std::array<T, Rows * Cols>::pointer;
    using const_pointer = typename std::array<T, Rows * Cols>::const_pointer;
    using reference = typename std::array<T, Rows * Cols>::reference;
    using const_reference = typename std::array<T, Rows * Cols>::const_reference;
    using iterator = typename std::array<T, Rows * Cols>::iterator;
    using const_iterator = typename std::array<T, Rows * Cols>::const_iterator;
    using size_type = typename std::array<T, Rows * Cols>::size_type;
    using difference_type = typename std::array<T, Rows * Cols>::difference_type;
    using reverse_iterator = typename std::array<T, Rows * Cols>::reverse_iterator;
    using const_reverse_iterator = typename std::array<T, Rows * Cols>::const_reverse_iterator;

  public:
    template <typename... U> constexpr Matrix(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == (Rows * Cols), "Wrong number of elements");
        static_assert((std::is_same_v<U, T> && ...), "All arguments must be exactly T");
    }

    constexpr Matrix(const std::array<T, Rows * Cols> &arr) : elements(arr)
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

    [[nodiscard]] constexpr inline bool operator==(const Matrix &v) const noexcept
        requires std::equality_comparable<T>
    {
        return elements == v.elements;
    }

    [[nodiscard]] constexpr inline bool operator!=(const Matrix &v) const noexcept
        requires std::equality_comparable<T>
    {
        return !(*this == v);
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return Rows * Cols;
    }

    std::array<T, Rows * Cols> elements{};
};

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols> &left,
                                                               const Matrix<T, Rows, Cols> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right.elements, std::plus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols> &left,
                                                               const T &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right, std::plus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator+(const T &left,
                                                               const Matrix<T, Rows, Cols> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left, right.elements, std::plus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols> &left,
                                                               const Matrix<T, Rows, Cols> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right.elements, std::minus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols> &left,
                                                               const T &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right, std::minus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator-(const T &left,
                                                               const Matrix<T, Rows, Cols> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left, right.elements, std::minus<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols> &v) noexcept
    requires Negatable<T>
{
    using namespace std;
    auto out = v;
    transform(begin(v), end(v), begin(out), std::negate<T>{});
    return out;
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator*(const Matrix<T, Rows, Cols> &left, T right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left.elements, right, std::multiplies<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator*(T left, const Matrix<T, Rows, Cols> &right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left, right.elements, std::multiplies<T>{});
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Matrix<T, Rows, Cols> operator/(const Matrix<T, Rows, Cols> &numerator,
                                                               T denominator) noexcept
    requires Divisable<T, T>
{
    return Elementwise(numerator.elements, denominator, std::divides<T>{});
}
