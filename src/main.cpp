#include <algorithm>
#include <cmath>
#include <iostream>

#include "image.hpp"
#include "types.hpp"

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

void DrawTriangle(Image &image, const Triangle &t, RGB color)
{
    BoundingBox box = TriangleBoundingBox(t);

    PointF32 top_left = {(f32)box.top_left.x, (f32)box.top_left.y};

    f32 inv_area_ABC = 1.0f / SignedTriangleArea(t);
    f32 area_PBC = SignedTriangleArea(top_left, t.p2, t.p3);
    f32 area_PAC = SignedTriangleArea(t.p1, top_left, t.p3);

    f32 u = area_PBC * inv_area_ABC;
    f32 v = area_PAC * inv_area_ABC;

    f32 uxscale = 0.5f * (t.p2.y - t.p3.y) * inv_area_ABC;
    f32 vxscale = 0.5f * (t.p3.y - t.p1.y) * inv_area_ABC;

    f32 uyscale = 0.5f * (t.p3.x - t.p2.x) * inv_area_ABC;
    f32 vyscale = 0.5f * (t.p1.x - t.p3.x) * inv_area_ABC;

    f32 u2 = (u - uxscale * (f32)box.top_left.x - uyscale * (f32)box.top_left.y);
    f32 v2 = (v - vxscale * (f32)box.top_left.x - vyscale * (f32)box.top_left.y);

    for (u32 y = 0; y < 400; ++y)
    {
        for (u32 x = 0; x < 400; ++x)
        {
            f32 a = u2 + (f32)y * uyscale + uxscale * (f32)x;
            f32 b = v2 + (f32)y * vyscale + vxscale * (f32)x;
            f32 c = 1.0f - a - b;

            if ((0 <= a) && (a <= 1.0f) && (0 <= b) && (b <= 1.0f) && (0 <= c) && (c <= 1.0f))
            {
                image[x, y] = color;
            }
        }
    }
}

inline void WritePixel(const RGB &color)
{
    std::cout << (int)color.r << ' ' << (int)color.g << ' ' << (int)color.b << '\n';
}

void WriteImage(const Image &image)
{
    using std::ranges::for_each;
    std::cout << "P3\n";
    std::cout << image.width() << ' ' << image.height() << "\n255\n";
    for_each(image, WritePixel);
}

int main()
{
    Triangle t1{

        {200, 100},
        {100, 300},
        {300, 300},
    };

    Triangle t2{

        {300, 100},
        {100, 300},
        {300, 300},
    };

    Image image(400, 400);
    DrawTriangle(image, t1, {255, 0, 0});
    DrawTriangle(image, t2, {0, 255, 0});
    WriteImage(image);
}
