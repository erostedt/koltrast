#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>
#include <vector>

#include "camera.hpp"
#include "counting_iterator.hpp"
#include "matrix.hpp"
#include "rasterizer.hpp"
#include "types.hpp"

using Mat4x4f = Matrix<f32, 4, 4>;
using Vec4f = Vector<f32, 4>;
using Vec3f = Vector<f32, 3>;
using Vec2f = Vector<f32, 2>;

template <typename T> struct BoundingBox
{
    Vector<T, 2> top_left;
    Vector<T, 2> bottom_right;
};

[[nodiscard]] constexpr inline size_t floor_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::floor(x));
}

[[nodiscard]] constexpr inline size_t ceil_to_size(f32 x) noexcept
{
    return static_cast<size_t>(std::ceil(x));
}

[[nodiscard]] constexpr inline bool is_empty(const BoundingBox<f32> &bb) noexcept
{
    return bb.top_left.x() >= bb.bottom_right.x() || bb.top_left.y() >= bb.bottom_right.y();
}

[[nodiscard]] constexpr inline BoundingBox<f32> bounding_box(const Vec4f &p1, const Vec4f &p2, const Vec4f &p3) noexcept
{
    using namespace std;
    f32 left = min(min(p1.x(), p2.x()), p3.x());
    f32 top = min(min(p1.y(), p2.y()), p3.y());
    f32 right = max(max(p1.x(), p2.x()), p3.x());
    f32 bottom = max(max(p1.y(), p2.y()), p3.y());
    return {{left, top}, {right, bottom}};
}

template <typename T> [[nodiscard]] constexpr inline BoundingBox<f32> bounding_box(const Image<T> &image) noexcept
{
    return {{0.0f, 0.0f}, {(f32)image.width(), (f32)image.height()}};
}

[[nodiscard]] constexpr inline BoundingBox<f32> intersect(const BoundingBox<f32> &bb1,
                                                          const BoundingBox<f32> &bb2) noexcept
{
    using namespace std;
    f32 l = max(bb1.top_left.x(), bb2.top_left.x());
    f32 r = min(bb1.bottom_right.x(), bb2.bottom_right.x());

    f32 t = max(bb1.top_left.y(), bb2.top_left.y());
    f32 b = min(bb1.bottom_right.y(), bb2.bottom_right.y());
    return {{l, t}, {r, b}};
}

[[nodiscard]] constexpr inline f32 edge(Vec4f p1, Vec4f p2, Vec4f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

[[nodiscard]] constexpr inline f32 edge(Vec2f p1, Vec2f p2, Vec2f p3) noexcept
{
    return (p2.y() - p1.y()) * (p3.x() - p1.x()) - (p2.x() - p1.x()) * (p3.y() - p1.y());
}

[[nodiscard]] constexpr inline bool out_of_z_bounds(const Vec4f &point) noexcept
{
    return point.z() < 0.0f || point.z() > 1.0f;
}

[[nodiscard]] constexpr inline BoundingBox<size_t> iteration_domain(const BoundingBox<f32> &bb) noexcept
{
    using namespace std;
    return {
        {
            floor_to_size(max(bb.top_left.x(), 0.0f)),
            floor_to_size(max(bb.top_left.y(), 0.0f)),
        },
        {
            ceil_to_size(bb.bottom_right.x()),
            ceil_to_size(bb.bottom_right.y()),
        },
    };
}

constexpr inline void rasterize_triangle(size_t triangle_index, const Vec4f &p1, const Vec4f &p2, const Vec4f &p3,
                                         const BoundingBox<f32> &bounds, DepthBuffer &depth_buffer,
                                         IndexBuffer &index_buffer) noexcept
{
    // Backface
    f32 tol = 1e-6f;
    f32 area = edge(p1, p2, p3);
    if (area <= tol)
    {
        return;
    }

    if (out_of_z_bounds(p1) && out_of_z_bounds(p2) && out_of_z_bounds(p3))
    {
        return;
    }

    const auto box = bounding_box(p1, p2, p3);
    const auto intersection = intersect(box, bounds);

    if (is_empty(intersection))
    {
        return;
    }

    const auto domain = iteration_domain(intersection);

    f32 left = (f32)domain.top_left.x() + 0.5f;
    f32 top = (f32)domain.top_left.y() + 0.5f;

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

    for (size_t y = domain.top_left.y(); y < domain.bottom_right.y(); ++y)
    {
        Vector<f32, 3> cw = weights;
        for (size_t x = domain.top_left.x(); x < domain.bottom_right.x(); ++x)
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

constexpr inline void _rasterize_triangles(const std::vector<Vec4f> &screen_vertices, const BoundingBox<f32> &bounds,
                                           DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    for (size_t i = 0; i < screen_vertices.size(); i += 3)
    {
        const auto v1 = screen_vertices[i + 0];
        const auto v2 = screen_vertices[i + 1];
        const auto v3 = screen_vertices[i + 2];
        rasterize_triangle(i / 3, v1, v2, v3, bounds, depth_buffer, index_buffer);
    }
}

constexpr inline void _rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                                           const BoundingBox<f32> &bounds, DepthBuffer &depth_buffer,
                                           IndexBuffer &index_buffer) noexcept
{
    for (size_t i = 0; i < faces.size(); ++i)
    {
        const auto &face = faces[i];
        const auto a = screen_vertices[face.vertex_indices[0]];
        const auto b = screen_vertices[face.vertex_indices[1]];
        const auto c = screen_vertices[face.vertex_indices[2]];
        rasterize_triangle(i, a, b, c, bounds, depth_buffer, index_buffer);
    }
}

constexpr inline RGB<f32> sample_bilinear(const Vec2f &uv, const Texture &texture) noexcept
{
    f32 x = std::clamp(uv.x() * (f32)(texture.width() - 1), 0.0f, (f32)(texture.width() - 1));
    f32 y = std::clamp(uv.y() * (f32)(texture.height() - 1), 0.0f, (f32)(texture.height() - 1));

    size_t fx = floor_to_size(x);
    size_t cx = ceil_to_size(x);
    size_t fy = floor_to_size(y);
    size_t cy = ceil_to_size(y);

    RGB<f32> a = texture[fx, fy];
    RGB<f32> b = texture[cx, fy];
    RGB<f32> c = texture[fx, cy];
    RGB<f32> d = texture[cx, cy];

    f32 s = x - (f32)fx;
    f32 t = y - (f32)fy;

    f32 red = (1.0f - s) * (1.0f - t) * a.r + s * (1.0f - t) * b.r + (1.0f - s) * t * c.r + s * t * d.r;
    f32 green = (1.0f - s) * (1.0f - t) * a.g + s * (1.0f - t) * b.g + (1.0f - s) * t * c.g + s * t * d.g;
    f32 blue = (1.0f - s) * (1.0f - t) * a.b + s * (1.0f - t) * b.b + (1.0f - s) * t * c.b + s * t * d.b;
    return {red, green, blue};
}

constexpr inline RGB<f32> sample_nearest_neighbor(const Vec2f &uv, const Texture &texture) noexcept
{
    size_t tx = floor_to_size(uv.x() * (f32)(texture.width() - 1));
    size_t ty = floor_to_size(uv.y() * (f32)(texture.height() - 1));
    return texture[tx, ty];
}

constexpr inline Vec2f interpolate_uv(const Vec2f &at, const Vec4f &v1, const Vec4f &v2, const Vec4f &v3,
                                      const Vec2f &uv1, const Vec2f &uv2, const Vec2f &uv3) noexcept
{
    Vec2f sv1 = {v1.x(), v1.y()};
    Vec2f sv2 = {v2.x(), v2.y()};
    Vec2f sv3 = {v3.x(), v3.y()};

    f32 w = 1.0f / edge(sv1, sv2, sv3);

    f32 ba = edge(sv2, sv3, at) * w;
    f32 bb = edge(sv3, sv1, at) * w;
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

    return uv_interp / invw_interp;
}

template <size_t Rows, size_t Cols>
constexpr inline Matrix<BoundingBox<f32>, Rows, Cols> make_grid(const Resolution &resolution) noexcept
{
    using namespace std;
    Matrix<BoundingBox<f32>, Rows, Cols> grid;
    for (size_t y = 0; y < Rows; ++y)
    {
        for (size_t x = 0; x < Cols; ++x)
        {
            f32 sx = (f32)x * (f32)resolution.width / (f32)Cols;
            f32 sy = (f32)y * (f32)resolution.height / (f32)Rows;
            f32 ex = (f32)(x + 1) * (f32)resolution.width / (f32)Cols;
            f32 ey = (f32)(y + 1) * (f32)resolution.height / (f32)Rows;
            grid[x, y] = {{sx, sy}, {ex, ey}};
        }
    }

    return grid;
}

/// ------------------------------------------ Public API ------------------------------------------

void reset_depth_buffer(DepthBuffer &buffer) noexcept
{
    using namespace std;
    fill(begin(buffer), end(buffer), numeric_limits<f32>::infinity());
}

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height)
{
    using namespace std;
    DepthBuffer buffer(width, height);
    reset_depth_buffer(buffer);
    return buffer;
}

void reset_index_buffer(IndexBuffer &buffer) noexcept
{
    using namespace std;
    fill(begin(buffer), end(buffer), numeric_limits<size_t>::max());
}

[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height)
{
    using namespace std;
    IndexBuffer buffer(width, height);
    reset_index_buffer(buffer);
    return buffer;
}

void rasterize_triangles(const std::vector<Vec4f> &screen_vertices, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept
{
    CHECK(screen_vertices.size() % 3 == 0);
    CHECK(depth_buffer.width() == index_buffer.width());
    CHECK(depth_buffer.height() == index_buffer.height());
    _rasterize_triangles(screen_vertices, bounding_box(depth_buffer), depth_buffer, index_buffer);
}

void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                         DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    CHECK(depth_buffer.width() == index_buffer.width());
    CHECK(depth_buffer.height() == index_buffer.height());

    using namespace std;
    const auto grid = make_grid<4, 4>({depth_buffer.width(), depth_buffer.height()});

    for_each(execution::par_unseq, begin(grid), end(grid), [&](const BoundingBox<f32> &bounds) {
        _rasterize_triangles(faces, screen_vertices, bounds, depth_buffer, index_buffer);
    });
}

void draw_triangles(ColorImage &image, const std::vector<RGB<u8>> &colors, const IndexBuffer &index_buffer) noexcept
{
    using namespace std;
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();

        size_t index = index_buffer[x, y];
        if (index != numeric_limits<size_t>::max())
        {
            image[x, y] = colors[index];
        }
    });
}

void draw_triangles(ColorImage &image, const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                    const std::vector<Vec2f> &texture_coordinates, const Texture &texture,
                    const IndexBuffer &index_buffer) noexcept
{
    using namespace std;
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();
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
            const Vec2f pixel_center = {(f32)x + 0.5f, (f32)y + 0.5f};
            const auto uv = interpolate_uv(pixel_center, v1, v2, v3, uv1, uv2, uv3);
            // image[x, y] = convert(sample_nearest_neighbor(uv, texture));
            image[x, y] = convert(sample_bilinear(uv, texture));
        }
    });
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
