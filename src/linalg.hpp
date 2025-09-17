#include "matrix.hpp"
#include "vector.hpp"

template <typename T, size_t Rows1, size_t Cols1, size_t Cols2>
[[nodiscard]] constexpr inline Matrix<T, Rows1, Cols2> operator*(const Matrix<T, Rows1, Cols1> &left,
                                                                 const Matrix<T, Cols1, Cols2> &right) noexcept
    requires Addable<T, T> && SelfAddable<T, T> && Multiplicable<T, T>
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
    requires Addable<T, T> && SelfAddable<T, T> && Multiplicable<T, T>
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
