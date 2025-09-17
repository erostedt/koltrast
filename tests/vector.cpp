#include "vector.hpp"

#include "types.hpp"
#include "utest.h"

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
