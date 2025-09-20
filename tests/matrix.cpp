#include "matrix.hpp"

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

UTEST(matrix, mat_pluseq_mat)
{
    Mat2 l = {
        {1, 2},
        {3, 4},
    };

    const Mat2 r = {
        {5, 6},
        {7, 8},
    };

    l += r;
    const Mat2 expected = {
        {6, 8},
        {10, 12},
    };
    ASSERT_EQ(expected, l);
}

UTEST(matrix, mat_pluseq_scalar)
{
    Mat2 l = {
        {1, 2},
        {3, 4},
    };

    const int scalar = 5;
    l += scalar;
    const Mat2 expected = {
        {6, 7},
        {8, 9},
    };
    ASSERT_EQ(expected, l);
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

UTEST(matrix, mat_minuseq_mat)
{
    Mat2 l = {
        {1, 2},
        {3, 4},
    };
    const Mat2 r = {
        {-1, 5},
        {6, 7},
    };
    l -= r;
    const Mat2 expected = {
        {2, -3},
        {-3, -3},
    };
    ASSERT_EQ(expected, l);
}

UTEST(matrix, mat_minuseq_scalar)
{
    Mat2 l = {
        {1, 2},
        {3, 6},
    };
    const int scalar = 5;
    l -= scalar;
    const Mat2 expected = {
        {-4, -3},
        {-2, 1},
    };
    ASSERT_EQ(expected, l);
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

UTEST(matrix, mat_timeseq_scalar)
{
    Mat2 l = {
        {1, 2},
        {3, 4},
    };
    const int scalar = 5;
    l *= scalar;
    const Mat2 expected = {
        {5, 10},
        {15, 20},
    };
    ASSERT_EQ(expected, l);
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

UTEST(matrix, mat_divedeeq_by_scalar)
{
    Mat2 l = {
        {5, 10},
        {-15, 20},
    };
    const int scalar = 5;
    l /= scalar;
    const Mat2 expected = {
        {1, 2},
        {-3, 4},
    };
    ASSERT_EQ(expected, l);
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
    const auto res = map(expected_eye - Mat::identity(), [](auto a) { return std::abs(a); });
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
    const auto res = map(expected_eye - Mat::identity(), [](auto a) { return std::abs(a); });
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

UTEST(matrix, transpose)
{
    using namespace std;
    const Matrix<i32, 3, 2> mat = {
        {1, 2},
        {3, 4},
        {5, 6},
    };
    const Matrix<i32, 2, 3> actual = mat.transposed();
    const Matrix<i32, 2, 3> expected = {{1, 3, 5}, {2, 4, 6}};
    ASSERT_EQ(expected, actual);
}

UTEST(vector2, accessors)
{
    const Vector<int, 2> u = {2, 1};
    ASSERT_EQ(2, u.x());
    ASSERT_EQ(1, u.y());
}

UTEST(vector3, accessors)
{
    const Vector<int, 3> u = {3, 2, 1};
    ASSERT_EQ(3, u.x());
    ASSERT_EQ(2, u.y());
    ASSERT_EQ(1, u.z());
}

UTEST(vector4, accessors)
{
    const Vector<int, 4> u = {4, 3, 2, 1};
    ASSERT_EQ(4, u.x());
    ASSERT_EQ(3, u.y());
    ASSERT_EQ(2, u.z());
    ASSERT_EQ(1, u.w());
}

UTEST(vector, vec_plus_vec)
{
    const Vector<int, 2> u = {1, 2};
    const Vector<int, 2> v = {-3, 4};
    const auto actual = u + v;
    const Vector<int, 2> expected = {-2, 6};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, vec_plus_scalar)
{
    const Vector<int, 2> u = {1, 2};
    const int scalar = 5;
    const auto actual = u + scalar;
    const Vector<int, 2> expected = {6, 7};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, scalar_plus_vec)
{
    const int scalar = 5;
    const Vector<int, 2> v = {1, 2};
    const auto actual = scalar + v;
    const Vector<int, 2> expected = {6, 7};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, vec_minus_vec)
{
    const Vector<int, 2> u = {1, 2};
    const Vector<int, 2> v = {-3, 4};
    const auto actual = u - v;
    const Vector<int, 2> expected = {4, -2};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, vec_minus_scalar)
{
    const Vector<int, 2> u = {1, 5};
    const int scalar = 5;
    const auto actual = u - scalar;
    const Vector<int, 2> expected = {-4, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, scalar_minus_vec)
{
    const int scalar = 5;
    const Vector<int, 2> u = {1, 5};
    const auto actual = scalar - u;
    const Vector<int, 2> expected = {4, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, negate_vec)
{
    const Vector<int, 2> u = {1, 5};
    const auto actual = -u;
    const Vector<int, 2> expected = {-1, -5};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, vec_times_scalar)
{
    const Vector<int, 2> u = {1, 5};
    const int scalar = 5;
    const auto actual = u * scalar;
    const Vector<int, 2> expected = {5, 25};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, scalar_times_vec)
{
    const int scalar = 5;
    const Vector<int, 2> u = {1, 5};
    const auto actual = scalar * u;
    const Vector<int, 2> expected = {5, 25};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, vec_diveded_by_scalar)
{
    const int scalar = 5;
    const Vector<int, 2> u = {25, 50};
    const auto actual = u / scalar;
    const Vector<int, 2> expected = {5, 10};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, equal)
{
    const Vector<int, 2> u = {1, 2};
    ASSERT_EQ(u, u);
}

UTEST(vector, not_equal)
{
    const Vector<int, 2> u = {0, 2};
    const Vector<int, 2> v = {1, 2};
    ASSERT_NE(u, v);
}

UTEST(vector, dot)
{
    const Vector<int, 3> u = {1, 2, 3};
    const Vector<int, 3> v = {-5, 3, 1};
    ASSERT_EQ(4, u.dot(v));
}

UTEST(vector, crossxy)
{
    const Vector<int, 3> u = {1, 0, 0};
    const Vector<int, 3> v = {0, 1, 0};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {0, 0, 1};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, crossyx)
{
    const Vector<int, 3> u = {0, 1, 0};
    const Vector<int, 3> v = {1, 0, 0};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {0, 0, -1};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, crossxz)
{
    const Vector<int, 3> u = {1, 0, 0};
    const Vector<int, 3> v = {0, 0, 1};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {0, -1, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, crosszx)
{
    const Vector<int, 3> u = {0, 0, 1};
    const Vector<int, 3> v = {1, 0, 0};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {0, 1, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, crossyz)
{
    const Vector<int, 3> u = {0, 1, 0};
    const Vector<int, 3> v = {0, 0, 1};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {1, 0, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, crosszy)
{
    const Vector<int, 3> u = {0, 0, 1};
    const Vector<int, 3> v = {0, 1, 0};
    const auto actual = u.cross(v);
    const Vector<int, 3> expected = {-1, 0, 0};
    ASSERT_EQ(expected, actual);
}

UTEST(vector, cross)
{
    const Vector<f32, 3> u = {1.0f, 2.0f, 3.0f};
    const Vector<f32, 3> v = {4.0f, 5.0f, 6.0f};
    const auto actual = u.cross(v);
    const Vector<f32, 3> expected = {-3.0f, 6.0f, -3.0f};
    ASSERT_NEAR(0.0f, (actual - expected).length(), 1e-5f);
}

UTEST(vector, size)
{
    const Vector<f32, 3> u = {1.0f, 2.0f, 3.0f};
    const auto actual = u.size();
    const size_t expected = 3;
    ASSERT_EQ(expected, actual);
}

UTEST(vector, squared_length)
{
    const Vector<f32, 3> u = {1.0f, 2.0f, 3.0f};
    const auto actual = u.squared_length();
    const f32 expected = 14.0f;
    ASSERT_NEAR(expected, actual, 1e-5f);
}

UTEST(vector, length)
{
    const Vector<f32, 3> u = {1.0f, 2.0f, 2.0f};
    const auto actual = u.length();
    const f32 expected = 3.0f;
    ASSERT_NEAR(expected, actual, 1e-5f);
}

UTEST(vector, squared_distance)
{
    const Vector<f32, 3> u = {1.0f, 1.0f, 1.0f};
    const Vector<f32, 3> v = {2.0f, 3.0f, 3.0f};
    const f32 actual = u.squared_distance(v);
    const f32 expected = 9.0f;
    ASSERT_NEAR(expected, actual, 1e-5f);
}

UTEST(vector, distance)
{
    const Vector<f32, 3> u = {1.0f, 1.0f, 1.0f};
    const Vector<f32, 3> v = {2.0f, 3.0f, 3.0f};
    const f32 actual = u.distance(v);
    const f32 expected = 3.0f;
    ASSERT_NEAR(expected, actual, 1e-5f);
}
