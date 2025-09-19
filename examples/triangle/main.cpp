#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include "camera.hpp"
#include "image.hpp"
#include "matrix.hpp"
#include "transform.hpp"
#include "types.hpp"

using Vec3f = Vector<f32, 3>;
using Vec2u = Vector<u32, 2>;
using ColorImage = Image<RGB>;
using DepthBuffer = Image<f32>;
using IndexBuffer = Image<size_t>;

struct Triangle
{
    Vec3f p1;
    Vec3f p2;
    Vec3f p3;
};

struct BoundingBox
{
    Vec2u top_left;
    Vec2u bottom_right;
};

inline u32 floor_to_u32(f32 x)
{
    return static_cast<u32>(std::floor(x));
}

inline u32 ceil_to_u32(f32 x)
{
    return static_cast<u32>(std::ceil(x));
}

BoundingBox clamped_triangle_bounding_box(const Triangle &t, const BoundingBox &bounds)
{
    assert(bounds.top_left.x() <= bounds.bottom_right.x());
    assert(bounds.top_left.y() <= bounds.bottom_right.y());

    using namespace std;
    f32 left = min(min(t.p1.x(), t.p2.x()), t.p3.x());
    u32 clamped_left = max(floor_to_u32(max(left, 0.0f)), bounds.top_left.x());

    f32 top = min(min(t.p1.y(), t.p2.y()), t.p3.y());
    u32 clamped_top = max(floor_to_u32(max(top, 0.0f)), bounds.top_left.y());

    f32 right = max(max(t.p1.x(), t.p2.x()), t.p3.x());
    u32 clamped_right = min(ceil_to_u32(max(right, 0.0f)), bounds.bottom_right.x());

    f32 bottom = max(max(t.p1.y(), t.p2.y()), t.p3.y());
    u32 clamped_bottom = min(ceil_to_u32(max(bottom, 0.0f)), bounds.bottom_right.y());
    return {{clamped_left, clamped_top}, {clamped_right, clamped_bottom}};
}

inline void write_pixel(const RGB &color)
{
    std::cout << (int)color.r << ' ' << (int)color.g << ' ' << (int)color.b << '\n';
}

void write_image(const ColorImage &image)
{
    using std::ranges::for_each;
    std::cout << "P3\n";
    std::cout << image.width() << ' ' << image.height() << "\n255\n";
    for_each(image, write_pixel);
}

inline f32 Edge(Vec3f p1, Vec3f p2, Vec3f p3)
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

inline f32 Edge(const Triangle &t)
{
    return Edge(t.p1, t.p2, t.p3);
}

inline DepthBuffer create_depth_buffer(u32 width, u32 height)
{
    using namespace std;
    DepthBuffer buffer(width, height);
    fill(begin(buffer), end(buffer), numeric_limits<f32>::infinity());
    return buffer;
}

inline IndexBuffer create_index_buffer(u32 width, u32 height)
{
    using namespace std;
    IndexBuffer buffer(width, height);
    fill(begin(buffer), end(buffer), numeric_limits<size_t>::max());
    return buffer;
}

inline bool out_of_bounds(const Vec3f &point, const BoundingBox &bounds)
{
    // TODO: CHECK IF CORRECT INEQUALITIES
    return point.x() < static_cast<f32>(bounds.top_left.x()) ||
           point.x() >= static_cast<f32>(bounds.bottom_right.x()) ||
           point.y() < static_cast<f32>(bounds.top_left.y()) ||
           point.y() >= static_cast<f32>(bounds.bottom_right.y()) || point.z() < 0.0f || point.z() > 1.0f;
}

void _rasterize_triangle(size_t triangle_index, const Triangle &t, DepthBuffer &depth_buffer, IndexBuffer &index_buffer)
{
    // Backface
    if (Edge(t) < 0.0f)
    {
        std::cout << "Backface" << std::endl;
        return;
    }

    BoundingBox bounds = {{0u, 0u}, {depth_buffer.width(), depth_buffer.height()}};

    if (out_of_bounds(t.p1, bounds) && out_of_bounds(t.p2, bounds) && out_of_bounds(t.p3, bounds))
    {
        std::cout << "OOB" << std::endl;
        std::cout << "P1: (" << t.p1.x() << " " << t.p1.y() << " " << t.p1.z() << ")" << std::endl;
        std::cout << "P2: (" << t.p2.x() << " " << t.p2.y() << " " << t.p2.z() << ")" << std::endl;
        std::cout << "P3: (" << t.p3.x() << " " << t.p3.y() << " " << t.p3.z() << ")" << std::endl;
        return;
    }

    BoundingBox box = clamped_triangle_bounding_box(t, bounds);

    f32 left = (f32)box.top_left.x() + 0.5f;
    f32 top = (f32)box.top_left.y() + 0.5f;

    f32 w = 1.0f / Edge(t.p1, t.p2, t.p3);
    f32 w0 = Edge(t.p2, t.p3, {left, top, 0.0f});
    f32 w1 = Edge(t.p3, t.p1, {left, top, 0.0f});
    f32 w2 = Edge(t.p1, t.p2, {left, top, 0.0f});

    f32 w0_dx = t.p3.y() - t.p2.y();
    f32 w0_dy = t.p2.x() - t.p3.x();

    f32 w1_dx = t.p1.y() - t.p3.y();
    f32 w1_dy = t.p3.x() - t.p1.x();

    f32 w2_dx = t.p2.y() - t.p1.y();
    f32 w2_dy = t.p1.x() - t.p2.x();

    for (u32 y = box.top_left.y(); y < box.bottom_right.y(); ++y)
    {
        f32 w0i = w0;
        f32 w1i = w1;
        f32 w2i = w2;

        for (u32 x = box.top_left.x(); x < box.bottom_right.x(); ++x)
        {
            if (w0i >= 0 && w1i >= 0 && w2i >= 0)
            {

                f32 l0 = w0i * w;
                f32 l1 = w1i * w;
                f32 l2 = w2i * w;

                f32 z = l0 * t.p1.z() + l1 * t.p2.z() + l2 * t.p3.z();
                if (z >= 0.0f && z <= 1.0f && z < depth_buffer[x, y])
                {
                    depth_buffer[x, y] = z;
                    index_buffer[x, y] = triangle_index;
                }
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

void rasterize_triangles(const std::vector<Triangle> &triangles, DepthBuffer &depth_buffer, IndexBuffer &index_buffer)
{
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        _rasterize_triangle(i, triangles[i], depth_buffer, index_buffer);
    }
}

void draw_triangles(ColorImage &image, const std::vector<RGB> &colors, const IndexBuffer &index_buffer)
{
    for (u32 y = 0; y < index_buffer.height(); ++y)
    {
        for (u32 x = 0; x < index_buffer.width(); ++x)
        {
            size_t index = index_buffer[x, y];
            if (index != std::numeric_limits<size_t>::max())
            {
                image[x, y] = colors[index];
            }
        }
    }
}

int main()
{
    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    const Vector<f32, 3> a{-1.0f, -1.0f, -1.0f};
    const Vector<f32, 3> b{2.0f, -1.0f, -3.0f};
    const Vector<f32, 3> c{0.0f, 1.0f, -2.0f};

    const Vector<f32, 3> d{1.0f, -1.0f, -1.0f};
    const Vector<f32, 3> e{0.0f, 1.0f, -2.0f};
    const Vector<f32, 3> f{-2.0f, -1.0f, -3.0f};

    const auto as = project_to_screen(a, view, proj, camera.resolution);
    const auto bs = project_to_screen(b, view, proj, camera.resolution);
    const auto cs = project_to_screen(c, view, proj, camera.resolution);

    const auto ds = project_to_screen(d, view, proj, camera.resolution);
    const auto es = project_to_screen(e, view, proj, camera.resolution);
    const auto fs = project_to_screen(f, view, proj, camera.resolution);

    Triangle t1{as, bs, cs};
    Triangle t2{ds, es, fs};

    ColorImage image((u32)camera.resolution.width, (u32)camera.resolution.height);
    auto depth_buffer = create_depth_buffer((u32)camera.resolution.width, (u32)camera.resolution.height);
    auto index_buffer = create_index_buffer((u32)camera.resolution.width, (u32)camera.resolution.height);
    rasterize_triangles({t1, t2}, depth_buffer, index_buffer);
    draw_triangles(image, {{255, 0, 0}, {0, 255, 0}}, index_buffer);
    write_image(image);
}
