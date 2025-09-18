#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "camera.hpp"
#include "image.hpp"
#include "matrix.hpp"
#include "transform.hpp"
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

    f32 right = max(max(t.p1.x, t.p2.x), t.p3.x);
    f32 bottom = max(max(t.p1.y, t.p2.y), t.p3.y);
    return {{FloorToU32(left), FloorToU32(top)}, {CeilToU32(right), CeilToU32(bottom)}};
}

BoundingBox ClampedTriangleBoundingBox(const Triangle &t, const BoundingBox &bounds)
{
    using namespace std;
    f32 left = min(min(min(t.p1.x, t.p2.x), t.p3.x), (f32)bounds.top_left.x);
    f32 top = min(min(min(t.p1.y, t.p2.y), t.p3.y), (f32)bounds.top_left.y);

    f32 right = max(max(max(t.p1.x, t.p2.x), t.p3.x), (f32)bounds.bottom_right.x);
    f32 bottom = max(max(max(t.p1.y, t.p2.y), t.p3.y), (f32)bounds.bottom_right.y);
    return {{FloorToU32(left), FloorToU32(top)}, {CeilToU32(right), CeilToU32(bottom)}};
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

inline f32 Edge(PointF32 p1, PointF32 p2, PointF32 p3)
{
    return (p2.y - p1.y) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.y - p1.y);
}

inline f32 Edge(const Triangle &t)
{
    return Edge(t.p1, t.p2, t.p3);
}

void DrawTriangle(Image &image, const Triangle &t, RGB color)
{
    assert(Edge(t) >= 0.0f && "Triangle must be Screen space CCW");
    BoundingBox box = ClampedTriangleBoundingBox(t, {{0, 0}, {image.width(), image.height()}});

    f32 left = (f32)box.top_left.x + 0.5f;
    f32 top = (f32)box.top_left.y + 0.5f;

    f32 w0 = Edge(t.p2, t.p3, {left, top});
    f32 w1 = Edge(t.p3, t.p1, {left, top});
    f32 w2 = Edge(t.p1, t.p2, {left, top});

    f32 w0_dx = t.p3.y - t.p2.y;
    f32 w0_dy = t.p2.x - t.p3.x;

    f32 w1_dx = t.p1.y - t.p3.y;
    f32 w1_dy = t.p3.x - t.p1.x;

    f32 w2_dx = t.p2.y - t.p1.y;
    f32 w2_dy = t.p1.x - t.p2.x;

    for (u32 y = box.top_left.y; y < box.bottom_right.y; ++y)
    {
        f32 w0i = w0;
        f32 w1i = w1;
        f32 w2i = w2;

        for (u32 x = box.top_left.x; x < box.bottom_right.x; ++x)
        {
            if (w0i >= 0 && w1i >= 0 && w2i >= 0)
            {
                image[x, y] = color;
            }
            w0i += w0_dx;
            w1i += w1_dx;
            w2i += w2_dx;
        }
        w0 += w0_dy;
        w1 += w1_dy;
        w2 += w2_dy;
    }
}

int main()
{
    const Vector<f32, 3> a{-1.0f, -1.0f, -2.0f};
    const Vector<f32, 3> b{1.0f, -1.0f, -2.0f};
    const Vector<f32, 3> c{0.0f, 1.0f, -2.0f};

    Camera<f32> camera = {{1280, 720}, 90, 1.0f, 10.0f};
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    const auto a1 = project_to_screen(a, view, proj, camera.resolution);
    const auto b1 = project_to_screen(b, view, proj, camera.resolution);
    const auto c1 = project_to_screen(c, view, proj, camera.resolution);

    Triangle t1{
        {a1.x(), a1.y()},
        {b1.x(), b1.y()},
        {c1.x(), c1.y()},
    };

    Image image((u32)camera.resolution.width, (u32)camera.resolution.height);
    DrawTriangle(image, t1, {255, 0, 0});
    WriteImage(image);
}
