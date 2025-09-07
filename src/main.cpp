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

inline f32 Edge(PointF32 p1, PointF32 p2, f32 x, f32 y)
{
    return (x - p1.x) * (p2.y - p1.y) - (y - p1.y) * (p2.x - p1.x);
}

void DrawTriangle(Image &image, const Triangle &t, RGB color)
{
    BoundingBox box = TriangleBoundingBox(t);
    for (u32 y = box.top_left.y; y < box.bottom_right.y; ++y)
    {
        for (u32 x = box.top_left.x; x < box.bottom_right.x; ++x)
        {
            f32 px = (f32)x + 0.5f;
            f32 py = (f32)y + 0.5f;
            f32 w0 = Edge(t.p2, t.p3, (f32)px, (f32)py);
            f32 w1 = Edge(t.p3, t.p1, (f32)px, (f32)py);
            f32 w2 = Edge(t.p1, t.p2, (f32)px, (f32)py);

            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
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
