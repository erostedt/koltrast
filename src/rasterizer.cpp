#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>
#include <vector>

#include "camera.hpp"
#include "counting_iterator.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "rasterizer.hpp"
#include "types.hpp"

template <typename T> struct BoundingBox
{
    Vector<T, 2> top_left;
    Vector<T, 2> bottom_right;
};

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
    f32 area = edge_function(p1, p2, p3);
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

    const Vec4f p = {(f32)domain.top_left.x() + 0.5f, (f32)domain.top_left.y() + 0.5f, 0.0f, 0.0f};

    f32 w = 1.0f / area;

    Vec3f weights = {
        edge_function(p2, p3, p),
        edge_function(p3, p1, p),
        edge_function(p1, p2, p),
    };

    Vec3f dx = {
        p3.y() - p2.y(),
        p1.y() - p3.y(),
        p2.y() - p1.y(),
    };

    Vec3f dy = {
        p2.x() - p3.x(),
        p3.x() - p1.x(),
        p1.x() - p2.x(),
    };

    Vec3f zs = {p1.z(), p2.z(), p3.z()};

    for (size_t y = domain.top_left.y(); y < domain.bottom_right.y(); ++y)
    {
        Vec3f cw = weights;
        for (size_t x = domain.top_left.x(); x < domain.bottom_right.x(); ++x)
        {
            if (cw.x() >= 0 && cw.y() >= 0 && cw.z() >= 0)
            {
                Vec3f lambdas = cw * w;

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

    using namespace std;
    const auto grid = make_grid<4, 4>({depth_buffer.width(), depth_buffer.height()});
    for_each(execution::par_unseq, begin(grid), end(grid), [&](const BoundingBox<f32> &bounds) {
        _rasterize_triangles(screen_vertices, bounds, depth_buffer, index_buffer);
    });
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

void dump_ppm(const ColorImage &image, std::ostream &stream)
{
    using std::ranges::for_each;
    stream << "P3\n";
    stream << image.width() << ' ' << image.height() << "\n255\n";

    const auto write_pixel = [&stream](const auto &color) {
        const auto r = std::clamp(color.r * 255.0f, 0.0f, 255.0f);
        const auto g = std::clamp(color.g * 255.0f, 0.0f, 255.0f);
        const auto b = std::clamp(color.b * 255.0f, 0.0f, 255.0f);
        stream << (int)r << ' ' << (int)g << ' ' << (int)b << '\n';
    };

    for_each(image, write_pixel);
}
