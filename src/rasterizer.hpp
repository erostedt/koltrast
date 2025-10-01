#pragma once
#include "camera.hpp"
#include "obj.hpp"
#include "types.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <execution>
#include <limits>
#include <vector>

#include "camera.hpp"
#include "check.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "types.hpp"

template <std::floating_point T> using ColorImage = Image<RGB<T>>;
template <std::floating_point T> using DepthBuffer = Image<T>;
using IndexBuffer = Image<size_t>;

template <typename T> struct BoundingBox
{
    Vec2<T> top_left;
    Vec2<T> bottom_right;

    [[nodiscard]] constexpr inline bool is_empty() const noexcept
    {
        return top_left.x() >= bottom_right.x() || top_left.y() >= bottom_right.y();
    }
};

template <std::floating_point T>
[[nodiscard]] constexpr inline BoundingBox<T> bounding_box(const Vec4<T> &p1, const Vec4<T> &p2,
                                                           const Vec4<T> &p3) noexcept
{
    using namespace std;
    T left = min(min(p1.x(), p2.x()), p3.x());
    T top = min(min(p1.y(), p2.y()), p3.y());
    T right = max(max(p1.x(), p2.x()), p3.x());
    T bottom = max(max(p1.y(), p2.y()), p3.y());
    return {{left, top}, {right, bottom}};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline BoundingBox<T> intersect(const BoundingBox<T> &bb1, const BoundingBox<T> &bb2) noexcept
{
    using namespace std;
    T l = max(bb1.top_left.x(), bb2.top_left.x());
    T r = min(bb1.bottom_right.x(), bb2.bottom_right.x());

    T t = max(bb1.top_left.y(), bb2.top_left.y());
    T b = min(bb1.bottom_right.y(), bb2.bottom_right.y());
    return {{l, t}, {r, b}};
}

template <std::floating_point T> [[nodiscard]] constexpr inline bool out_of_z_bounds(const Vec4<T> &point) noexcept
{
    return point.z() < T{0} || point.z() > T{1};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline BoundingBox<size_t> iteration_domain(const BoundingBox<T> &bb) noexcept
{
    using namespace std;
    return {
        {
            floor_to_size(max(bb.top_left.x(), T{0})),
            floor_to_size(max(bb.top_left.y(), T{0})),
        },
        {
            ceil_to_size(bb.bottom_right.x()),
            ceil_to_size(bb.bottom_right.y()),
        },
    };
}

template <std::floating_point T> constexpr inline bool inside_triangle(const Vec3<T> &bary) noexcept
{
    return bary.x() >= 0 && bary.y() >= 0 && bary.z() >= 0;
}

template <std::floating_point T>
constexpr inline void rasterize_triangle(size_t triangle_index, const Vec4<T> &p1, const Vec4<T> &p2, const Vec4<T> &p3,
                                         const BoundingBox<T> &bounds, DepthBuffer<T> &depth_buffer,
                                         IndexBuffer &index_buffer) noexcept
{
    // Backface
    T tol = T(1e-6);
    T area = edge_function(p1, p2, p3);
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

    if (intersection.is_empty())
    {
        return;
    }

    const auto domain = iteration_domain(intersection);

    const Vec4<T> p = {(T)domain.top_left.x() + T{0.5}, (T)domain.top_left.y() + T{0.5}, 0.0f, 0.0f};

    T w = T{1.0} / area;

    Vec3<T> weights = {
        edge_function(p2, p3, p),
        edge_function(p3, p1, p),
        edge_function(p1, p2, p),
    };

    Vec3<T> dx = {
        p3.y() - p2.y(),
        p1.y() - p3.y(),
        p2.y() - p1.y(),
    };

    Vec3<T> dy = {
        p2.x() - p3.x(),
        p3.x() - p1.x(),
        p1.x() - p2.x(),
    };

    Vec3<T> zs = {p1.z(), p2.z(), p3.z()};

    for (size_t y = domain.top_left.y(); y < domain.bottom_right.y(); ++y)
    {
        Vec3<T> cw = weights;
        for (size_t x = domain.top_left.x(); x < domain.bottom_right.x(); ++x)
        {
            if (inside_triangle(cw))
            {
                Vec3<T> lambdas = cw * w;

                T z = lambdas.dot(zs);
                if (z >= T{0} && z <= T{1} && z < depth_buffer[x, y])
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

template <std::floating_point T>
constexpr inline void _rasterize_triangles(const std::vector<Vec4<T>> &screen_vertices, const BoundingBox<T> &bounds,
                                           DepthBuffer<T> &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    for (size_t i = 0; i < screen_vertices.size(); i += 3)
    {
        const auto v1 = screen_vertices[i + 0];
        const auto v2 = screen_vertices[i + 1];
        const auto v3 = screen_vertices[i + 2];
        rasterize_triangle(i / 3, v1, v2, v3, bounds, depth_buffer, index_buffer);
    }
}

template <std::floating_point T>
constexpr inline void _rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4<T>> &screen_vertices,
                                           const BoundingBox<T> &bounds, DepthBuffer<T> &depth_buffer,
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

template <std::floating_point T, size_t Rows, size_t Cols>
constexpr inline Matrix<BoundingBox<T>, Rows, Cols> make_grid(const Resolution &resolution) noexcept
{
    using namespace std;
    Matrix<BoundingBox<T>, Rows, Cols> grid;
    for (size_t y = 0; y < Rows; ++y)
    {
        for (size_t x = 0; x < Cols; ++x)
        {
            T sx = (T)x * (T)resolution.width / (T)Cols;
            T sy = (T)y * (T)resolution.height / (T)Rows;
            T ex = min((T)(x + 1) * (T)resolution.width / (T)Cols, resolution.width - 1);
            T ey = min((T)(y + 1) * (T)resolution.height / (T)Rows, resolution.height - 1);
            grid[x, y] = {{sx, sy}, {ex, ey}};
        }
    }

    return grid;
}

/// ------------------------------------------ Public API ------------------------------------------
template <std::floating_point T> constexpr inline void reset_depth_buffer(DepthBuffer<T> &buffer) noexcept
{
    using namespace std;
    fill(begin(buffer), end(buffer), numeric_limits<T>::infinity());
}
template <std::floating_point T> [[nodiscard]] inline DepthBuffer<T> create_depth_buffer(size_t width, size_t height)
{
    using namespace std;
    DepthBuffer<T> buffer(width, height);
    reset_depth_buffer(buffer);
    return buffer;
}
constexpr inline void reset_index_buffer(IndexBuffer &buffer) noexcept
{
    using namespace std;
    fill(begin(buffer), end(buffer), numeric_limits<size_t>::max());
}
[[nodiscard]] inline IndexBuffer create_index_buffer(size_t width, size_t height)
{
    using namespace std;
    IndexBuffer buffer(width, height);
    reset_index_buffer(buffer);
    return buffer;
}

template <std::floating_point T, size_t RowTiles = 4, size_t ColTiles = RowTiles>
    requires(RowTiles > 0) && (ColTiles > 0)
inline void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4<T>> &screen_vertices,
                                DepthBuffer<T> &depth_buffer, IndexBuffer &index_buffer) noexcept
{
    CHECK(depth_buffer.width() == index_buffer.width());
    CHECK(depth_buffer.height() == index_buffer.height());

    using namespace std;
    const auto grid = make_grid<T, RowTiles, ColTiles>({depth_buffer.width(), depth_buffer.height()});
    for_each(execution::par_unseq, begin(grid), end(grid), [&](const BoundingBox<T> &bounds) {
        _rasterize_triangles(faces, screen_vertices, bounds, depth_buffer, index_buffer);
    });
}

template <std::floating_point T> inline void dump_ppm(const ColorImage<T> &image, std::ostream &stream)
{
    using std::ranges::for_each;
    stream << "P3\n";
    stream << image.width() << ' ' << image.height() << "\n255\n";

    const auto write_pixel = [&stream](const auto &color) {
        const auto r = std::clamp(color.r * T{255}, T{0}, T{255});
        const auto g = std::clamp(color.g * T{255}, T{0}, T{255});
        const auto b = std::clamp(color.b * T{255}, T{0}, T{255});
        stream << (int)r << ' ' << (int)g << ' ' << (int)b << '\n';
    };

    for_each(image, write_pixel);
}
