#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>

#include "camera.hpp"
#include "matrix.hpp"
#include "rasterizer.hpp"
#include "types.hpp"

using Mat4x4f = Matrix<f32, 4, 4>;
using Vec4f = Vector<f32, 4>;
using Vec3f = Vector<f32, 3>;
using Vec2f = Vector<f32, 2>;

[[nodiscard]] constexpr inline size_t floor_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::floor(x));
}

[[nodiscard]] constexpr inline size_t ceil_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::ceil(x));
}

[[nodiscard]] constexpr inline BoundingBox clamped_triangle_bounding_box(const Vec4f &p1, const Vec4f &p2,
                                                                         const Vec4f &p3,
                                                                         const BoundingBox &bounds) noexcept
{
    CHECK(bounds.top_left.x() <= bounds.bottom_right.x());
    CHECK(bounds.top_left.y() <= bounds.bottom_right.y());

    using namespace std;
    f32 left = min(min(p1.x(), p2.x()), p3.x());
    size_t clamped_left = max(floor_to_size(max(left, 0.0f)), bounds.top_left.x());

    f32 top = min(min(p1.y(), p2.y()), p3.y());
    size_t clamped_top = max(floor_to_size(max(top, 0.0f)), bounds.top_left.y());

    f32 right = max(max(p1.x(), p2.x()), p3.x());
    size_t clamped_right = min(ceil_to_size(max(right, 0.0f)), bounds.bottom_right.x());

    f32 bottom = max(max(p1.y(), p2.y()), p3.y());
    size_t clamped_bottom = min(ceil_to_size(max(bottom, 0.0f)), bounds.bottom_right.y());
    return {{clamped_left, clamped_top}, {clamped_right, clamped_bottom}};
}

[[nodiscard]] constexpr inline f32 edge(Vec4f p1, Vec4f p2, Vec4f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

[[nodiscard]] constexpr inline f32 edge(Vec2f p1, Vec2f p2, Vec2f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

[[nodiscard]] constexpr inline bool out_of_bounds(const Vec4f &point, const BoundingBox &bounds) noexcept
{
    // TODO: CHECK IF CORRECT INEQUALITIES
    return point.x() < static_cast<f32>(bounds.top_left.x()) ||
           point.x() >= static_cast<f32>(bounds.bottom_right.x()) ||
           point.y() < static_cast<f32>(bounds.top_left.y()) ||
           point.y() >= static_cast<f32>(bounds.bottom_right.y()) || point.z() < 0.0f || point.z() > 1.0f;
}

inline void rasterize_triangle(size_t triangle_index, const Vec4f &p1, const Vec4f &p2, const Vec4f &p3,
                               DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    // Backface
    f32 area = edge(p1, p2, p3);
    if (area < 0.0f)
    {
        std::cerr << "Backface" << std::endl;
        return;
    }

    BoundingBox bounds = {{0u, 0u}, {depth_buffer.width(), depth_buffer.height()}};

    if (out_of_bounds(p1, bounds) && out_of_bounds(p2, bounds) && out_of_bounds(p3, bounds))
    {
        return;
    }

    BoundingBox box = clamped_triangle_bounding_box(p1, p2, p3, bounds);

    f32 left = (f32)box.top_left.x() + 0.5f;
    f32 top = (f32)box.top_left.y() + 0.5f;

    f32 w = 1.0f / area;

    Vector<f32, 3> weights = {
        edge(p2, p3, {left, top, 0.0f, 0.0f}),
        edge(p3, p1, {left, top, 0.0f, 0.0f}),
        edge(p1, p2, {left, top, 0.0f, 0.0f}),
    };

    Vector<f32, 3> dx = {
        p3.y() - p2.y(),
        p1.y() - p3.y(),
        p2.y() - p1.y(),
    };

    Vector<f32, 3> dy = {
        p2.x() - p3.x(),
        p3.x() - p1.x(),
        p1.x() - p2.x(),
    };

    Vector<f32, 3> zs = {p1.z(), p2.z(), p3.z()};

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

void rasterize_triangles(const std::vector<Vec4f> &screen_vertices, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept
{
    CHECK(screen_vertices.size() % 3 == 0);
    for (size_t i = 0; i < screen_vertices.size(); i += 3)
    {
        const auto v1 = screen_vertices[i + 0];
        const auto v2 = screen_vertices[i + 1];
        const auto v3 = screen_vertices[i + 2];
        rasterize_triangle(i / 3, v1, v2, v3, depth_buffer, index_buffer);
    }
}

void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                         DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    for (size_t i = 0; i < faces.size(); ++i)
    {
        const auto &face = faces[i];
        const auto a = screen_vertices[face.vertex_indices[0]];
        const auto b = screen_vertices[face.vertex_indices[1]];
        const auto c = screen_vertices[face.vertex_indices[2]];
        rasterize_triangle(i, a, b, c, depth_buffer, index_buffer);
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

inline RGB sample_texture_nearest_neighbor(const Image<RGB> &texture, size_t x, size_t y, const Vec4f &v1,
                                           const Vec4f &v2, const Vec4f &v3, const Vec2f &uv1, const Vec2f &uv2,
                                           const Vec2f &uv3)
{
    Vec2f pixel_center = {static_cast<f32>(x) + 0.5f, static_cast<f32>(y) + 0.5f};
    Vec2f sv1 = {v1.x(), v1.y()};
    Vec2f sv2 = {v2.x(), v2.y()};
    Vec2f sv3 = {v3.x(), v3.y()};

    f32 w = 1.0f / edge(sv1, sv2, sv3);

    f32 ba = edge(sv2, sv3, pixel_center) * w;
    f32 bb = edge(sv3, sv1, pixel_center) * w;
    f32 bc = 1.0f - ba - bb;

    f32 wa = v1.w();
    f32 wb = v2.w();
    f32 wc = v3.w();

    f32 iwa = 1.0f / wa;
    f32 iwb = 1.0f / wb;
    f32 iwc = 1.0f / wc;

    Vec2f uv1i = uv1 * iwa;
    Vec2f uv2i = uv2 * iwb;
    Vec2f uv3i = uv3 * iwc;

    Vec2f uv_interp = ba * uv1i + bb * uv2i + bc * uv3i;
    f32 invw_interp = ba * iwa + bb * iwb + bc * iwc;

    Vec2f uv = uv_interp / invw_interp;

    // Nearest neighbor
    size_t tx = floor_to_size(std::clamp(uv.x(), 0.0f, 0.999999f) * (f32)texture.width());
    size_t ty = floor_to_size(std::clamp(uv.y(), 0.0f, 0.999999f) * (f32)texture.height());
    return texture[tx, ty];
}

void draw_triangles(ColorImage &image, const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                    const std::vector<Vec2f> &texture_coordinates, const Image<RGB> &texture,
                    const IndexBuffer &index_buffer) noexcept
{
    for (size_t y = 0; y < index_buffer.height(); ++y)
    {
        for (size_t x = 0; x < index_buffer.width(); ++x)
        {
            size_t index = index_buffer[x, y];
            if (index != std::numeric_limits<size_t>::max())
            {
                const auto &face = faces[index];
                const auto v1 = screen_vertices[face.vertex_indices[0]];
                const auto v2 = screen_vertices[face.vertex_indices[1]];
                const auto v3 = screen_vertices[face.vertex_indices[2]];

                const auto uv1 = texture_coordinates[face.texture_indices[0]];
                const auto uv2 = texture_coordinates[face.texture_indices[1]];
                const auto uv3 = texture_coordinates[face.texture_indices[2]];
                image[x, y] = sample_texture_nearest_neighbor(texture, x, y, v1, v2, v3, uv1, uv2, uv3);
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
