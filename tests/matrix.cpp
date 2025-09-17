#include "matrix.hpp"

#include "utest.h"

using Mat2 = Matrix<int, 2, 2>;

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
