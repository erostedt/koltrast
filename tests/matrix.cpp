#include "matrix.hpp"
#include "linalg.hpp"

#include "utest.h"

using Mat2 = Matrix<int, 2, 2>;

UTEST(matrix, mat_index)
{
    const Matrix<int, 3, 2> m = {1, 2, 3, 4, 5, 6};
    ASSERT_EQ(m.index(0, 0), 0ul);
    ASSERT_EQ(m.index(0, 1), 1ul);
    ASSERT_EQ(m.index(1, 0), 2ul);
    ASSERT_EQ(m.index(1, 1), 3ul);
    ASSERT_EQ(m.index(2, 0), 4ul);
    ASSERT_EQ(m.index(2, 1), 5ul);
}

UTEST(matrix, mat_access_1d)
{
    const Matrix<int, 3, 2> m = {1, 2, 3, 4, 5, 6};
    ASSERT_EQ(m[0], 1);
    ASSERT_EQ(m[1], 2);
    ASSERT_EQ(m[2], 3);
    ASSERT_EQ(m[3], 4);
    ASSERT_EQ(m[4], 5);
    ASSERT_EQ(m[5], 6);
}

UTEST(matrix, mat_access_2d)
{
    const Matrix<int, 3, 2> m = {1, 2, 3, 4, 5, 6};
    const auto m00 = m[0, 0];
    const auto m01 = m[0, 1];
    const auto m10 = m[1, 0];
    const auto m11 = m[1, 1];
    const auto m20 = m[2, 0];
    const auto m21 = m[2, 1];
    ASSERT_EQ(m00, 1);
    ASSERT_EQ(m01, 2);
    ASSERT_EQ(m10, 3);
    ASSERT_EQ(m11, 4);
    ASSERT_EQ(m20, 5);
    ASSERT_EQ(m21, 6);
}

UTEST(matrix, mat_sizes)
{
    const Matrix<int, 3, 2> m = {1, 2, 3, 4, 5, 6};
    ASSERT_EQ(m.size(), 6ul);
    ASSERT_EQ(m.rows(), 3ul);
    ASSERT_EQ(m.cols(), 2ul);
}

UTEST(matrix, mat_plus_mat)
{
    const Mat2 l = {1, 2, 3, 4};
    const Mat2 r = {5, 6, 7, 8};

    const auto actual = l + r;
    const Mat2 expected = {6, 8, 10, 12};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_plus_scalar)
{
    const Mat2 l = {1, 2, 3, 4};
    const int scalar = 5;
    const auto actual = l + scalar;
    const Mat2 expected = {6, 7, 8, 9};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, scalar_plus_mat)
{
    const int scalar = 5;
    const Mat2 r = {1, 2, 3, 4};
    const auto actual = scalar + r;
    const Mat2 expected = {6, 7, 8, 9};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_minus_mat)
{
    const Mat2 l = {1, 2, 3, 4};
    const Mat2 r = {-1, 5, 6, 7};
    const auto actual = l - r;
    const Mat2 expected = {2, -3, -3, -3};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_minus_scalar)
{
    const Mat2 l = {1, 2, 3, 6};
    const int scalar = 5;
    const auto actual = l - scalar;
    const Mat2 expected = {-4, -3, -2, 1};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, scalar_minus_mat)
{
    const int scalar = 5;
    const Mat2 r = {1, 2, 3, 6};
    const auto actual = scalar - r;
    const Mat2 expected = {4, 3, 2, -1};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, negate_mat)
{
    const Mat2 m = {-1, 1, -2, 2};
    const auto actual = -m;
    const Mat2 expected = {1, -1, 2, -2};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_times_scalar)
{
    const Mat2 l = {1, 2, 3, 4};
    const int scalar = 5;
    const auto actual = l * scalar;
    const Mat2 expected = {5, 10, 15, 20};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, scalar_times_vec)
{
    const int scalar = 5;
    const Mat2 r = {1, 2, 3, 4};
    const auto actual = scalar * r;
    const Mat2 expected = {5, 10, 15, 20};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_diveded_by_scalar)
{
    const Mat2 l = {5, 10, -15, 20};
    const int scalar = 5;
    const auto actual = l / scalar;
    const Mat2 expected = {1, 2, -3, 4};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, equal)
{
    const Mat2 l = {5, 10, -15, 20};
    ASSERT_EQ(l, l);
}

UTEST(matrix, not_equal)
{
    const Mat2 l = {1, 2, 3, 4};
    const Mat2 r = {0, 2, 3, 4};
    ASSERT_NE(l, r);
}

UTEST(matrix, mat_mul_mat)
{
    Matrix<int, 2, 3> l = {1, 2, 3, 4, 5, 6};
    Matrix<int, 3, 2> r = {7, 8, 9, 10, 11, 12};
    const auto actual = l * r;
    const Matrix<int, 2, 2> expected = {58, 64, 139, 154};
    ASSERT_EQ(actual, expected);
}

UTEST(matrix, mat_mul_vec)
{
    Matrix<int, 2, 3> l = {1, 2, 3, 4, 5, 6};
    Vector<int, 3> r = {7, 8, 9};
    const auto actual = l * r;
    const Vector<int, 2> expected = {50, 122};
    ASSERT_EQ(actual, expected);
}
