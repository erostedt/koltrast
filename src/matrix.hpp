#pragma once

#include "elementwise.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <expected>
#include <functional>
#include <initializer_list>
#include <string>

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
    using Mat = Matrix<T, Rows, Cols>;

  public:
    constexpr Matrix() = default;

    template <typename... U> constexpr Matrix(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == (Rows * Cols), "Wrong number of elements");
        static_assert((std::is_convertible_v<U, T> && ...), "All arguments must be exactly T");
    }

    constexpr Matrix(const std::array<T, Rows * Cols> &arr) : elements(arr)
    {
    }

    Matrix(std::initializer_list<std::initializer_list<T>> mat) : elements{}
    {
        assert(std::size(mat) == Rows);
        size_t index = 0;
        for (const auto &row : mat)
        {
            assert(std::size(row) == Cols);
            for (const T &element : row)
            {
                elements[index] = element;
                ++index;
            }
        }
    }

    constexpr static Mat identity() noexcept
    {
        Mat eye;
        for (size_t i = 0; i < Rows; ++i)
        {
            eye[i, i] = T{1};
        }
        return eye;
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

    [[nodiscard]] static constexpr inline size_t rows() noexcept

    {
        return Rows;
    }

    [[nodiscard]] static constexpr inline size_t cols() noexcept
    {
        return Cols;
    }

    [[nodiscard]] static constexpr inline size_t size() noexcept
    {
        return Rows * Cols;
    }

    [[nodiscard]] constexpr inline size_t index(size_t row, size_t col) const noexcept
    {
        return row * Cols + col;
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t index) const noexcept
    {
        return elements[index];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t index) noexcept
    {
        return elements[index];
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t row, size_t col) const noexcept
    {
        return elements[index(row, col)];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t row, size_t col) noexcept
    {
        return elements[index(row, col)];
    }

    [[nodiscard]] constexpr std::expected<Mat, std::string> inverse() const noexcept
        requires std::floating_point<T> && (Rows == Cols)
    {
        Mat A = *this;
        Mat inv = Mat::identity();

        // Gauss–Jordan elimination
        for (size_t col = 0; col < Cols; ++col)
        {
            // Find pivot row
            size_t pivot = col;
            T max_val = std::abs(A[pivot, col]);
            for (size_t row = col + 1; row < Rows; ++row)
            {
                T val = std::abs(A[row, col]);
                if (val > max_val)
                {
                    pivot = row;
                    max_val = val;
                }
            }

            if (max_val < std::numeric_limits<T>::epsilon())
            {
                return std::unexpected("Matrix is singular and cannot be inverted");
            }

            // Swap rows in both matrices if needed
            if (pivot != col)
            {
                for (size_t j = 0; j < Cols; ++j)
                {
                    std::swap(A[col, j], A[pivot, j]);
                    std::swap(inv[col, j], inv[pivot, j]);
                }
            }

            // Normalize pivot row
            T diag = A[col, col];
            for (size_t j = 0; j < Cols; ++j)
            {
                A[col, j] /= diag;
                inv[col, j] /= diag;
            }

            // Eliminate other rows
            for (size_t row = 0; row < Rows; ++row)
            {
                if (row == col)
                {
                    continue;
                }
                T factor = A[row, col];
                for (size_t j = 0; j < Cols; ++j)
                {
                    A[row, j] -= factor * A[col, j];
                    inv[row, j] -= factor * inv[col, j];
                }
            }
        }

        return std::expected<Mat, std::string>(inv);
    }

    template <typename UnaryOperation> [[nodiscard]] constexpr inline auto elementwise(UnaryOperation op) const noexcept
    {
        using out_type = decltype(op(std::declval<T>()));
        return Matrix<out_type, Rows, Cols>(Elementwise(elements, op));
    }

    template <typename U, typename BinaryOperation>
    [[nodiscard]] constexpr inline auto elementwise(const Matrix<U, Rows, Cols> &other,
                                                    BinaryOperation op) const noexcept
    {
        using out_type = decltype(op(std::declval<T>(), std::declval<U>()));
        return Matrix<out_type, Rows, Cols>(Elementwise(elements, other, op));
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
