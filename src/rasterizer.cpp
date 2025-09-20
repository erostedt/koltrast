#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "camera.hpp"
#include "matrix.hpp"
#include "rasterizer.hpp"
#include "types.hpp"

using Vec3f = Vector<f32, 3>;

[[nodiscard]] constexpr inline size_t floor_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::floor(x));
}

[[nodiscard]] constexpr inline size_t ceil_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::ceil(x));
}

[[nodiscard]] constexpr inline BoundingBox clamped_triangle_bounding_box(const Triangle &t,
                                                                         const BoundingBox &bounds) noexcept
{
    CHECK(bounds.top_left.x() <= bounds.bottom_right.x());
    CHECK(bounds.top_left.y() <= bounds.bottom_right.y());

    using namespace std;
    f32 left = min(min(t.p1.x(), t.p2.x()), t.p3.x());
    size_t clamped_left = max(floor_to_size(max(left, 0.0f)), bounds.top_left.x());

    f32 top = min(min(t.p1.y(), t.p2.y()), t.p3.y());
    size_t clamped_top = max(floor_to_size(max(top, 0.0f)), bounds.top_left.y());

    f32 right = max(max(t.p1.x(), t.p2.x()), t.p3.x());
    size_t clamped_right = min(ceil_to_size(max(right, 0.0f)), bounds.bottom_right.x());

    f32 bottom = max(max(t.p1.y(), t.p2.y()), t.p3.y());
    size_t clamped_bottom = min(ceil_to_size(max(bottom, 0.0f)), bounds.bottom_right.y());
    return {{clamped_left, clamped_top}, {clamped_right, clamped_bottom}};
}

[[nodiscard]] constexpr inline f32 edge(Vec3f p1, Vec3f p2, Vec3f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

[[nodiscard]] constexpr inline f32 edge(const Triangle &t) noexcept
{
    return edge(t.p1, t.p2, t.p3);
}

[[nodiscard]] constexpr inline bool out_of_bounds(const Vec3f &point, const BoundingBox &bounds) noexcept
{
    // TODO: CHECK IF CORRECT INEQUALITIES
    return point.x() < static_cast<f32>(bounds.top_left.x()) ||
           point.x() >= static_cast<f32>(bounds.bottom_right.x()) ||
           point.y() < static_cast<f32>(bounds.top_left.y()) ||
           point.y() >= static_cast<f32>(bounds.bottom_right.y()) || point.z() < 0.0f || point.z() > 1.0f;
}

constexpr inline void rasterize_triangle(size_t triangle_index, const Triangle &t, DepthBuffer &depth_buffer,
                                         IndexBuffer &index_buffer) noexcept
{
    // Backface
    if (edge(t) < 0.0f)
    {
        std::cerr << "Backface" << std::endl;
        return;
    }

    BoundingBox bounds = {{0u, 0u}, {depth_buffer.width(), depth_buffer.height()}};

    if (out_of_bounds(t.p1, bounds) && out_of_bounds(t.p2, bounds) && out_of_bounds(t.p3, bounds))
    {
        std::cerr << "OOB" << std::endl;
        std::cerr << "P1: (" << t.p1.x() << " " << t.p1.y() << " " << t.p1.z() << ")" << std::endl;
        std::cerr << "P2: (" << t.p2.x() << " " << t.p2.y() << " " << t.p2.z() << ")" << std::endl;
        std::cerr << "P3: (" << t.p3.x() << " " << t.p3.y() << " " << t.p3.z() << ")" << std::endl;
        return;
    }

    BoundingBox box = clamped_triangle_bounding_box(t, bounds);

    f32 left = (f32)box.top_left.x() + 0.5f;
    f32 top = (f32)box.top_left.y() + 0.5f;

    f32 w = 1.0f / edge(t.p1, t.p2, t.p3);

    Vector<f32, 3> weights = {
        edge(t.p2, t.p3, {left, top, 0.0f}),
        edge(t.p3, t.p1, {left, top, 0.0f}),
        edge(t.p1, t.p2, {left, top, 0.0f}),
    };

    Vector<f32, 3> dx = {
        t.p3.y() - t.p2.y(),
        t.p1.y() - t.p3.y(),
        t.p2.y() - t.p1.y(),
    };

    Vector<f32, 3> dy = {
        t.p2.x() - t.p3.x(),
        t.p3.x() - t.p1.x(),
        t.p1.x() - t.p2.x(),
    };

    Vector<f32, 3> zs = {t.p1.z(), t.p2.z(), t.p3.z()};

    for (size_t y = box.top_left.y(); y < box.bottom_right.y(); ++y)
    {
        Vector<f32, 3> cw = weights;
        for (size_t x = box.top_left.x(); x < box.bottom_right.x(); ++x)
        {
            if (cw.x() >= 0 && cw.y() >= 0 && cw.z() >= 0)
            {
                Vector<f32, 3> lambdas = cw * w;

                f32 z = lambdas.dot(zs);
                if (z >= 0.0f && z <= 1.0f && z < depth_buffer[x, y])
                {
                    depth_buffer[x, y] = z;
                    index_buffer[x, y] = triangle_index;
                }
            }
            cw += dx;
        }
        weights += dy;
    }
}

/// ------------------------------------------ Public API ------------------------------------------

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height)
{
    using namespace std;
    DepthBuffer buffer(width, height);
    fill(begin(buffer), end(buffer), numeric_limits<f32>::infinity());
    return buffer;
}

[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height)
{
    using namespace std;
    IndexBuffer buffer(width, height);
    fill(begin(buffer), end(buffer), numeric_limits<size_t>::max());
    return buffer;
}

void rasterize_triangles(const std::vector<Triangle> &triangles, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept
{
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        rasterize_triangle(i, triangles[i], depth_buffer, index_buffer);
    }
}

void draw_triangles(ColorImage &image, const std::vector<RGB> &colors, const IndexBuffer &index_buffer) noexcept
{
    for (size_t y = 0; y < index_buffer.height(); ++y)
    {
        for (size_t x = 0; x < index_buffer.width(); ++x)
        {
            size_t index = index_buffer[x, y];
            if (index != std::numeric_limits<size_t>::max())
            {
                image[x, y] = colors[index];
            }
        }
    }
}

void dump_ppm(const ColorImage &image, std::ostream &stream)
{
    using std::ranges::for_each;
    stream << "P3\n";
    stream << image.width() << ' ' << image.height() << "\n255\n";

    const auto write_pixel = [&stream](const auto &color) {
        stream << (int)color.r << ' ' << (int)color.g << ' ' << (int)color.b << '\n';
    };

    for_each(image, write_pixel);
}
