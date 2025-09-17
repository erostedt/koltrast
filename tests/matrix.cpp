#include "matrix.hpp"
#include "linalg.hpp"

#include "types.hpp"
#include "utest.h"

using Mat2 = Matrix<int, 2, 2>;

UTEST(matrix, mat_index)
{
    const Matrix<int, 3, 2> m = {
        {1, 2},
        {3, 4},
        {5, 6},
    };

    ASSERT_EQ(0ul, m.index(0, 0));
    ASSERT_EQ(1ul, m.index(0, 1));
    ASSERT_EQ(2ul, m.index(1, 0));
    ASSERT_EQ(3ul, m.index(1, 1));
    ASSERT_EQ(4ul, m.index(2, 0));
    ASSERT_EQ(5ul, m.index(2, 1));
}

UTEST(matrix, mat_access_1d)
{
    const Matrix<int, 3, 2> m = {
        {1, 2},
        {3, 4},
        {5, 6},
    };

    ASSERT_EQ(1, m[0]);
    ASSERT_EQ(2, m[1]);
    ASSERT_EQ(3, m[2]);
    ASSERT_EQ(4, m[3]);
    ASSERT_EQ(5, m[4]);
    ASSERT_EQ(6, m[5]);
}

UTEST(matrix, mat_access_2d)
{
    const Matrix<int, 3, 2> m = {
        {1, 2},
        {3, 4},
        {5, 6},
    };

    const auto m00 = m[0, 0];
    const auto m01 = m[0, 1];
    const auto m10 = m[1, 0];
    const auto m11 = m[1, 1];
    const auto m20 = m[2, 0];
    const auto m21 = m[2, 1];
    ASSERT_EQ(1, m00);
    ASSERT_EQ(2, m01);
    ASSERT_EQ(3, m10);
    ASSERT_EQ(4, m11);
    ASSERT_EQ(5, m20);
    ASSERT_EQ(6, m21);
}

UTEST(matrix, mat_sizes)
{
    const Matrix<int, 3, 2> m = {
        {1, 2},
        {3, 4},
        {5, 6},
    };

    ASSERT_EQ(6ul, m.size());
    ASSERT_EQ(3ul, m.rows());
    ASSERT_EQ(2ul, m.cols());
}

UTEST(matrix, mat_plus_mat)
{
    const Mat2 l = {
        {1, 2},
        {3, 4},
    };

    const Mat2 r = {
        {5, 6},
        {7, 8},
    };

    const auto actual = l + r;
    const Mat2 expected = {
        {6, 8},
        {10, 12},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_plus_scalar)
{
    const Mat2 l = {
        {1, 2},
        {3, 4},
    };

    const int scalar = 5;
    const auto actual = l + scalar;
    const Mat2 expected = {
        {6, 7},
        {8, 9},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, scalar_plus_mat)
{
    const int scalar = 5;
    const Mat2 r = {
        {1, 2},
        {3, 4},
    };

    const auto actual = scalar + r;
    const Mat2 expected = {
        {6, 7},
        {8, 9},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_minus_mat)
{
    const Mat2 l = {
        {1, 2},
        {3, 4},
    };
    const Mat2 r = {
        {-1, 5},
        {6, 7},
    };
    const auto actual = l - r;
    const Mat2 expected = {
        {2, -3},
        {-3, -3},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_minus_scalar)
{
    const Mat2 l = {
        {1, 2},
        {3, 6},
    };
    const int scalar = 5;
    const auto actual = l - scalar;
    const Mat2 expected = {
        {-4, -3},
        {-2, 1},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, scalar_minus_mat)
{
    const int scalar = 5;
    const Mat2 r = {
        {1, 2},
        {3, 6},
    };
    const auto actual = scalar - r;
    const Mat2 expected = {
        {4, 3},
        {2, -1},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, negate_mat)
{
    const Mat2 m = {
        {-1, 1},
        {-2, 2},
    };
    const auto actual = -m;
    const Mat2 expected = {
        {1, -1},
        {2, -2},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_times_scalar)
{
    const Mat2 l = {
        {1, 2},
        {3, 4},
    };
    const int scalar = 5;
    const auto actual = l * scalar;
    const Mat2 expected = {
        {5, 10},
        {15, 20},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, scalar_times_vec)
{
    const int scalar = 5;
    const Mat2 r = {
        {1, 2},
        {3, 4},
    };
    const auto actual = scalar * r;
    const Mat2 expected = {
        {5, 10},
        {15, 20},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_diveded_by_scalar)
{
    const Mat2 l = {
        {5, 10},
        {-15, 20},
    };
    const int scalar = 5;
    const auto actual = l / scalar;
    const Mat2 expected = {
        {1, 2},
        {-3, 4},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, equal)
{
    const Mat2 l = {
        {5, 10},
        {-15, 20},
    };
    ASSERT_EQ(l, l);
}

UTEST(matrix, not_equal)
{
    const Mat2 l = {
        {1, 2},
        {3, 4},
    };
    const Mat2 r = {
        {0, 2},
        {3, 4},
    };
    ASSERT_NE(l, r);
}

UTEST(matrix, mat_mul_mat)
{
    Matrix<int, 2, 3> l = {
        {1, 2, 3},
        {4, 5, 6},
    };
    Matrix<int, 3, 2> r = {
        {7, 8},
        {9, 10},
        {11, 12},
    };
    const auto actual = l * r;
    const Matrix<int, 2, 2> expected = {
        {58, 64},
        {139, 154},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_mul_vec)
{
    Matrix<int, 2, 3> l = {
        {1, 2, 3},
        {4, 5, 6},
    };
    Matrix<int, 2, 3>::identity();
    Vector<int, 3> r = {7, 8, 9};
    const auto actual = l * r;
    const Vector<int, 2> expected = {50, 122};
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_identity)
{
    using namespace std;
    using Mat = Matrix<int, 4, 4>;
    const Mat actual = Mat::identity();
    const Mat expected = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
    };
    ASSERT_EQ(expected, actual);
}

UTEST(matrix, mat_inverse1)
{
    using namespace std;
    using Mat = Matrix<f64, 4, 4>;
    Mat mat = {
        {0.832921, 0.0, 0.553392, 0.0},
        {0.291613, 0.849893, -0.438913, 0.0},
        {-0.470323, 0.526956, 0.707894, 0.0},
        {-2.574104, 3.650642, 4.868381, 1.0},
    };
    const auto inv = mat.inverse();
    ASSERT_TRUE(inv.has_value());
    const auto expected_eye = mat * inv.value();
    const auto res = (expected_eye - Mat::identity()).elementwise([](auto a) { return std::abs(a); });
    const auto sad = reduce(begin(res), end(res), 0.0);

    ASSERT_NEAR(0.0, sad, 1e-5);
}

UTEST(matrix, mat_inverse2)
{
    using namespace std;
    using Mat = Matrix<f64, 4, 4>;
    const Mat mat = {
        {4.0, 7.0, 2.0, 3.0},
        {0.0, 5.0, 0.0, 1.0},
        {0.0, 0.0, 3.0, 0.0},
        {0.0, 0.0, 0.0, 2.0},
    };

    const auto inv = mat.inverse();
    ASSERT_TRUE(inv.has_value());
    const auto expected_eye = mat * inv.value();
    const auto res = (expected_eye - Mat::identity()).elementwise([](auto a) { return std::abs(a); });
    const auto sad = reduce(begin(res), end(res), 0.0);

    ASSERT_NEAR(0.0, sad, 1e-5);
}

UTEST(matrix, mat_inverse_singular_rows_dependant)
{
    using namespace std;
    using Mat = Matrix<f64, 4, 4>;
    const Mat mat = {
        {1.0, 2.0, 3.0, 4.0},
        {2.0, 4.0, 6.0, 8.0},
        {0.0, 1.0, 0.0, 1.0},
        {0.0, 0.0, 1.0, 1.0},
    };
    const auto inv = mat.inverse();
    ASSERT_FALSE(inv.has_value());
}

UTEST(matrix, mat_inverse_singular_cols_dependant)
{
    using namespace std;
    using Mat = Matrix<f64, 4, 4>;
    const Mat mat = {
        {1.0, 2.0, 3.0, 4.0},
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 0.0, 0.0, 0.0},
        {5.0, 6.0, 7.0, 8.0},
    };
    const auto inv = mat.inverse();
    ASSERT_FALSE(inv.has_value());
}
