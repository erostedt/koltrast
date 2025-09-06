#include <cmath>
#include <iostream>

using u32 = uint32_t;
using f32 = float;

template <typename T> struct Point
{
    T x;
    T y;
};

using PointF32 = Point<f32>;
using PointU32 = Point<u32>;

struct Triangle
{
    PointF32 p1;
    PointF32 p2;
    PointF32 p3;
};

struct BoundingBox
{
    PointU32 top_left;
    PointU32 bottom_right;
};

inline f32 SignedTriangleArea(PointF32 p1, PointF32 p2, PointF32 p3)
{
    f32 det = p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y);
    return 0.5f * det;
}

inline f32 SignedTriangleArea(const Triangle &t)
{
    return SignedTriangleArea(t.p1, t.p2, t.p3);
}

inline u32 FloorToU32(f32 x)
{
    return static_cast<u32>(std::floor(x));
}

inline u32 CeilToU32(f32 x)
{
    return static_cast<u32>(std::ceil(x));
}

BoundingBox TriangleBoundingBox(const Triangle &t)
{
    using namespace std;
    f32 left = min(min(t.p1.x, t.p2.x), t.p3.x);
    f32 top = min(min(t.p1.y, t.p2.y), t.p3.y);

    f32 right = max(max(t.p1.y, t.p2.y), t.p3.y);
    f32 bottom = max(max(t.p1.x, t.p2.x), t.p3.x);
    return {{FloorToU32(left), FloorToU32(top)}, {CeilToU32(right), CeilToU32(bottom)}};
}

int main()
{
    Triangle triangle{

        {200, 100},
        {100, 300},
        {300, 300},
    };

    BoundingBox box = TriangleBoundingBox(triangle);

    PointF32 top_left = {(f32)box.top_left.x, (f32)box.top_left.y};

    f32 inv_area_ABC = 1.0f / SignedTriangleArea(triangle);
    f32 area_PBC = SignedTriangleArea(top_left, triangle.p2, triangle.p3);
    f32 area_PAC = SignedTriangleArea(triangle.p1, top_left, triangle.p3);

    f32 u = area_PBC * inv_area_ABC;
    f32 v = area_PAC * inv_area_ABC;

    f32 uxscale = 0.5f * (triangle.p2.y - triangle.p3.y) * inv_area_ABC;
    f32 vxscale = 0.5f * (triangle.p3.y - triangle.p1.y) * inv_area_ABC;

    f32 uyscale = 0.5f * (triangle.p3.x - triangle.p2.x) * inv_area_ABC;
    f32 vyscale = 0.5f * (triangle.p1.x - triangle.p3.x) * inv_area_ABC;

    std::cout << "P3\n";
    std::cout << 400 << ' ' << 400 << "\n255\n";
    for (u32 y = 0; y < 400; ++y)
    {
        f32 dy = (f32)(y - box.top_left.y);
        for (u32 x = 0; x < 400; ++x)
        {
            f32 dx = (f32)(x - box.top_left.x);

            f32 a = u + dy * uyscale + dx * uxscale;
            f32 b = v + dy * vyscale + dx * vxscale;
            f32 c = 1.0f - a - b;

            if ((0 <= a) && (a <= 1.0f) && (0 <= b) && (b <= 1.0f) && (0 <= c) && (c <= 1.0f))
            {
                std::cout << "255 0 0\n";
            }
            else
            {
                std::cout << "0 0 0\n";
            }
        }
    }
}
