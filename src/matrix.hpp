#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <expected>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <string>

template <typename T>
concept Addable = requires(T x) {
    { x + x } -> std::same_as<T>;
};

template <typename T>
concept SelfAddable = requires(T x) {
    { x += x };
};

template <typename T>
concept Subtractable = requires(T x) {
    { x - x } -> std::same_as<T>;
};

template <typename T>
concept Negatable = requires(T x) {
    { -x } -> std::same_as<T>;
};

template <typename T>
concept Multiplicable = requires(T x) {
    { x * x } -> std::same_as<T>;
};

template <typename T>
concept Divisable = requires(T a) {
    { a / a } -> std::same_as<T>;
};

template <typename T>
concept FixedSize = requires {
    { std::remove_cvref_t<T>::size() } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept ElementwiseOperable = std::ranges::input_range<T> && FixedSize<T>;

template <typename T>
concept ElementwiseAddable = ElementwiseOperable<T> && Addable<typename T::value_type>;

template <typename T>
concept ElementwiseSubtractable = ElementwiseOperable<T> && Subtractable<typename T::value_type>;

template <typename T>
concept ElementwiseMultiplicable = ElementwiseOperable<T> && Multiplicable<typename T::value_type>;

template <typename T>
concept ElementwiseDivisable = ElementwiseOperable<T> && Divisable<typename T::value_type>;

template <typename T>
concept ElementwiseNegatable = ElementwiseOperable<T> && Negatable<typename T::value_type>;

template <typename T>
concept ElementwiseEquatable = ElementwiseOperable<T> && std::equality_comparable<typename T::value_type>;

template <ElementwiseOperable Container, typename UnaryOperation>
constexpr Container Elementwise(const Container &left, UnaryOperation op)
{
    using namespace std;
    Container out;
    transform(begin(left), end(left), begin(out), op);
    return out;
}

template <ElementwiseOperable Container, typename BinaryOperation>
constexpr Container Elementwise(const Container &left, const Container &right, BinaryOperation op)
{
    using namespace std;
    Container out;
    transform(begin(left), end(left), begin(right), begin(out), op);
    return out;
}

template <ElementwiseOperable Container, typename BinaryOperation>
constexpr Container Elementwise(const Container &left, const typename Container::value_type &right, BinaryOperation op)
{
    using namespace std;
    using T = typename Container::value_type;
    Container out;
    transform(begin(left), end(left), begin(out), [&right, &op](const T &element) { return op(element, right); });
    return out;
}

template <ElementwiseOperable Container, typename BinaryOperation>
constexpr Container Elementwise(const typename Container::value_type &left, const Container &right, BinaryOperation op)
{
    using namespace std;
    using T = typename Container::value_type;
    Container out;
    transform(begin(right), end(right), begin(out), [&left, &op](const T &element) { return op(left, element); });
    return out;
}

template <ElementwiseAddable Container>
[[nodiscard]] constexpr inline Container operator+(const Container &left, const Container &right) noexcept
{
    return Elementwise(left, right, std::plus<typename Container::value_type>{});
}

template <ElementwiseAddable Container>
[[nodiscard]] constexpr inline Container operator+(const Container &left,
                                                   const typename Container::value_type &right) noexcept
{
    return Elementwise(left, right, std::plus<typename Container::value_type>{});
}

template <ElementwiseAddable Container>
[[nodiscard]] constexpr inline Container operator+(const typename Container::value_type &left,
                                                   const Container &right) noexcept
{
    return Elementwise(left, right, std::plus<typename Container::value_type>{});
}

template <ElementwiseSubtractable Container>
[[nodiscard]] constexpr inline Container operator-(const Container &left, const Container &right) noexcept
{
    return Elementwise(left, right, std::minus<typename Container::value_type>{});
}

template <ElementwiseSubtractable Container>
[[nodiscard]] constexpr inline Container operator-(const Container &left,
                                                   const typename Container::value_type &right) noexcept
{
    return Elementwise(left, right, std::minus<typename Container::value_type>{});
}

template <ElementwiseSubtractable Container>
[[nodiscard]] constexpr inline Container operator-(const typename Container::value_type &left,
                                                   const Container &right) noexcept
{
    return Elementwise(left, right, std::minus<typename Container::value_type>{});
}

template <ElementwiseNegatable Container>
[[nodiscard]] constexpr inline Container operator-(const Container &v) noexcept
{
    using namespace std;
    auto out = v;
    transform(begin(v), end(v), begin(out), std::negate<typename Container::value_type>{});
    return out;
}

template <ElementwiseMultiplicable Container>
[[nodiscard]] constexpr inline Container operator*(const Container &left,
                                                   const typename Container::value_type &right) noexcept
{
    return Elementwise(left, right, std::multiplies<typename Container::value_type>{});
}

template <ElementwiseMultiplicable Container>
[[nodiscard]] constexpr inline Container operator*(const typename Container::value_type &left,
                                                   const Container &right) noexcept
{
    return Elementwise(left, right, std::multiplies<typename Container::value_type>{});
}

template <ElementwiseDivisable Container>
[[nodiscard]] constexpr inline Container operator/(const Container &numerator,
                                                   const typename Container::value_type &denominator) noexcept
{
    return Elementwise(numerator, denominator, std::divides<typename Container::value_type>{});
}

template <ElementwiseEquatable Container>
[[nodiscard]] constexpr inline bool operator==(const Container &left, const Container &right) noexcept
{
    using namespace std;
    return equal(begin(left), end(left), begin(right));
}

template <ElementwiseEquatable Container>
[[nodiscard]] constexpr inline bool operator!=(const Container &left, const Container &right) noexcept
{
    return !(left == right);
}

template <typename T, size_t N>
    requires(N > 0)
class Array
{
  public:
    using value_type = typename std::array<T, N>::value_type;
    using pointer = typename std::array<T, N>::pointer;
    using const_pointer = typename std::array<T, N>::const_pointer;
    using reference = typename std::array<T, N>::reference;
    using const_reference = typename std::array<T, N>::const_reference;
    using iterator = typename std::array<T, N>::iterator;
    using const_iterator = typename std::array<T, N>::const_iterator;
    using size_type = typename std::array<T, N>::size_type;
    using difference_type = typename std::array<T, N>::difference_type;
    using reverse_iterator = typename std::array<T, N>::reverse_iterator;
    using const_reverse_iterator = typename std::array<T, N>::const_reverse_iterator;

    constexpr Array() = default;
    template <typename... U> constexpr Array(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == N, "Wrong number of elements");
        static_assert((std::is_convertible_v<U, T> && ...), "All arguments must be exactly T");
    }
    constexpr Array(const std::array<T, N> &arr) : elements{arr}
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

    [[nodiscard]] static consteval inline size_t size() noexcept
    {
        return N;
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t index) const noexcept
    {
        return elements[index];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t index) noexcept
    {
        return elements[index];
    }

  protected:
    std::array<T, N> elements{};
};

template <typename Matrix>
concept Square = requires(Matrix m) {
    { m.rows() == m.cols() };
};

template <typename T, size_t Rows, size_t Cols>
    requires std::is_signed_v<T>
class Matrix : public Array<T, Rows * Cols>
{
  public:
    using Array<T, Rows * Cols>::Array;
    using Array<T, Rows * Cols>::elements;
    using Array<T, Rows * Cols>::operator[];

    Matrix(std::initializer_list<std::initializer_list<T>> mat)
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

    [[nodiscard]] constexpr static Matrix identity() noexcept
    {
        Matrix eye;
        for (size_t i = 0; i < Rows; ++i)
        {
            eye[i, i] = T{1};
        }
        return eye;
    }

    [[nodiscard]] static constexpr inline size_t rows() noexcept
    {
        return Rows;
    }

    [[nodiscard]] static constexpr inline size_t cols() noexcept
    {
        return Cols;
    }

    [[nodiscard]] constexpr inline size_t index(size_t row, size_t col) const noexcept
    {
        return row * Cols + col;
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t row, size_t col) const noexcept
    {
        return elements[index(row, col)];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t row, size_t col) noexcept
    {
        return elements[index(row, col)];
    }

    [[nodiscard]] constexpr std::expected<Matrix, std::string> inverse() const noexcept
        requires std::floating_point<T> && (Rows == Cols)
    {
        Matrix A = *this;
        Matrix inv = Matrix::identity();

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

        return std::expected<Matrix, std::string>(inv);
    }
};

template <typename T, size_t Rows>
    requires std::is_signed_v<T>
class Vector : public Array<T, Rows>
{
    using Array<T, Rows>::Array;
    using Array<T, Rows>::elements;

  public:
    [[nodiscard]] constexpr inline T dot(const Vector &v) const noexcept
        requires Multiplicable<T> && Addable<T>
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

    [[nodiscard]] constexpr inline Vector cross(const Vector &v) const noexcept
        requires(Rows == 3) && Multiplicable<T> && Subtractable<T>
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
};

template <typename T, size_t Rows1, size_t Cols1, size_t Cols2>
[[nodiscard]] constexpr inline Matrix<T, Rows1, Cols2> operator*(const Matrix<T, Rows1, Cols1> &left,
                                                                 const Matrix<T, Cols1, Cols2> &right) noexcept
    requires Addable<T> && SelfAddable<T> && Multiplicable<T>
{
    Matrix<T, Rows1, Cols2> result;
    for (size_t row = 0; row < Rows1; ++row)
    {
        for (size_t col = 0; col < Cols2; ++col)
        {
            T accum{};
            for (size_t k = 0; k < Cols1; ++k)
            {
                accum += left[row, k] * right[k, col];
            }
            result[row, col] = accum;
        }
    }

    return result;
}

template <typename T, size_t Rows, size_t Cols>
[[nodiscard]] constexpr inline Vector<T, Rows> operator*(const Matrix<T, Rows, Cols> &left,
                                                         const Vector<T, Cols> &right) noexcept
    requires Addable<T> && SelfAddable<T> && Multiplicable<T>
{
    Vector<T, Rows> result;
    for (size_t row = 0; row < Rows; ++row)
    {
        T accum{};
        for (size_t col = 0; col < Cols; ++col)
        {
            accum += left[row, col] * right[col];
        }
        result[row] = accum;
    }

    return result;
}
