#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <limits>
#include <vector>

#include "check.hpp"
#include "counting_iterator.hpp"
#include "depth_buffer.hpp"
#include "image.hpp"
#include "math.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "shader.hpp"
#include "types.hpp"

using IndexBuffer = Buffer3<size_t>;

template <std::floating_point T> inline std::vector<Vec2<T>> make_aa_grid(size_t rows, size_t cols)
{
    CHECK(rows > 0);
    CHECK(cols > 0);

    std::vector<Vec2<T>> offsets;
    offsets.reserve(rows * cols);

    T x_delta = T{1} / T(cols);
    T y_delta = T{1} / T(rows);
    for (size_t y = 0; y < rows; ++y)
    {
        for (size_t x = 0; x < cols; ++x)
        {
            T fx = (T(x) + T(0.5)) * x_delta - T(0.5);
            T fy = (T(y) + T(0.5)) * y_delta - T(0.5);
            offsets.push_back({fx, fy});
        }
    }
    return offsets;
}

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
[[nodiscard]] constexpr inline BoundingBox<T> bounding_box(const Vec3<T> &p1, const Vec3<T> &p2,
                                                           const Vec3<T> &p3) noexcept
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

template <std::floating_point T> [[nodiscard]] constexpr inline bool out_of_z_bounds(const Vec3<T> &point) noexcept
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
            T ex = min((T)(x + 1) * (T)resolution.width / (T)Cols, (T)resolution.width - T{1});
            T ey = min((T)(y + 1) * (T)resolution.height / (T)Rows, (T)resolution.height - T{1});
            grid[x, y] = {{sx, sy}, {ex, ey}};
        }
    }

    return grid;
}

template <std::floating_point T>
constexpr inline Vec3<T> barycentric(const Vec2<T> &at, const Vec3<T> &v1, const Vec3<T> &v2,
                                     const Vec3<T> &v3) noexcept
{
    const Vec3<T> p = {at.x(), at.y(), T{0}};
    T w = T{1} / edge_function(v1, v2, v3);
    T a = edge_function(v2, v3, p) * w;
    T b = edge_function(v3, v1, p) * w;
    return {a, b, T{1} - a - b};
}

template <std::floating_point T>
constexpr inline void rasterize_pixel(size_t triangle_index, const Vec3<T> &weights, const T w, const Vec3<T> &dx,
                                      const Vec3<T> &dy, const Vec3<T> &zs, const std::vector<Vec2<T>> &sample_offsets,
                                      std::span<T> depth_cell, std::span<size_t> index_cell) noexcept
{
    for (size_t i = 0; i < sample_offsets.size(); ++i)
    {
        const Vec2<T> offset = sample_offsets[i];
        const Vec3<T> offsetted = weights + dx * offset.x() + dy * offset.y();
        if (inside_triangle(offsetted))
        {
            Vec3<T> lambdas = offsetted * w;

            T z = lambdas.dot(zs);
            if (z >= T{0} && z <= T{1} && z < depth_cell[i])
            {
                depth_cell[i] = z;
                index_cell[i] = triangle_index;
            }
        }
    }
}

template <std::floating_point T>
constexpr inline void rasterize_triangle(size_t triangle_index, const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &p3,
                                         const BoundingBox<T> &bounds, const std::vector<Vec2<T>> &sample_offsets,
                                         DepthBuffer<T> &depth_buffer, IndexBuffer &index_buffer) noexcept
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

    const Vec3<T> p = {(T)domain.top_left.x() + T{0.5}, (T)domain.top_left.y() + T{0.5}, 0.0f};

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
            rasterize_pixel(triangle_index, cw, w, dx, dy, zs, sample_offsets, depth_buffer.slice(x, y),
                            index_buffer.slice(x, y));
            cw += dx;
        }
        weights += dy;
    }
}

template <std::floating_point T, size_t TileRows = 8, size_t TileCols = TileRows>
inline void rasterize_triangles(const std::vector<Vec3<T>> &screen_coordinates,
                                const std::vector<Vec2<T>> &sample_offsets, DepthBuffer<T> &depth_buffer,
                                IndexBuffer &index_buffer) noexcept
{
    CHECK(screen_coordinates.size() % 3 == 0);
    CHECK(depth_buffer.width() == index_buffer.width());
    CHECK(depth_buffer.height() == index_buffer.height());
    CHECK(depth_buffer.depth() == index_buffer.depth());

    using namespace std;
    const auto grid = make_grid<T, TileRows, TileCols>({depth_buffer.width(), depth_buffer.height()});
    for_each(execution::par_unseq, begin(grid), end(grid), [&](const BoundingBox<T> &bounds) {
        for (size_t i = 0; i < screen_coordinates.size() / 3; ++i)
        {
            const Vec3<T> &a = screen_coordinates[3 * i + 0];
            const Vec3<T> &b = screen_coordinates[3 * i + 1];
            const Vec3<T> &c = screen_coordinates[3 * i + 2];
            rasterize_triangle(i, a, b, c, bounds, sample_offsets, depth_buffer, index_buffer);
        }
    });
}

template <std::floating_point T, FragmentShader<T> FragmentShader>
inline void render_vertices(const std::vector<Vec3<T>> &screen_coordinates,
                            const std::vector<OutputVertex<T>> &vertex_buffer, const FragmentShader &fragment_shader,
                            const std::vector<Vec2<T>> &sample_offsets, const IndexBuffer &index_buffer,
                            ColorImage<T> &linear_image) noexcept
{
    CHECK(index_buffer.depth() == sample_offsets.size());
    using namespace std;
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();
        std::span<const size_t> indicies = index_buffer.slice(x, y);
        if (all_of(begin(indicies), end(indicies), [](size_t idx) { return idx == numeric_limits<size_t>::max(); }))
        {
            return;
        }

        const T sample_contribution = T{1} / T(indicies.size());

        size_t coverage = 0;
        RGB<T> color = BLACK<T>;
        for (size_t d = 0; d < indicies.size(); ++d)
        {
            size_t index = indicies[d];
            if (index != numeric_limits<size_t>::max())
            {
                const Vec3<T> sv1 = screen_coordinates[3 * index + 0];
                const Vec3<T> sv2 = screen_coordinates[3 * index + 1];
                const Vec3<T> sv3 = screen_coordinates[3 * index + 2];

                const OutputVertex<T> v1 = vertex_buffer[3 * index + 0];
                const OutputVertex<T> v2 = vertex_buffer[3 * index + 1];
                const OutputVertex<T> v3 = vertex_buffer[3 * index + 2];

                const Vec2<T> offset = sample_offsets[d];

                T px = (T)x + T(0.5) + offset.x();
                T py = (T)y + T(0.5) + offset.y();

                const Vec3<T> bary = barycentric({px, py}, sv1, sv2, sv3);
                const Vec3<T> world_position = bary.x() * v1.position + bary.y() * v2.position + bary.z() * v3.position;
                const Vec3<T> world_normal = bary.x() * v1.normal + bary.y() * v2.normal + bary.z() * v3.normal;
                const Vec2<T> uv = texture_uv(bary, v1.clip_position.w(), v2.clip_position.w(), v3.clip_position.w(),
                                              v1.texture_coordinates, v2.texture_coordinates, v3.texture_coordinates);

                color = color + sample_contribution * fragment_shader({world_position, world_normal, uv});
                ++coverage;
            }
        }
        T coverage_fraction = (T)coverage * sample_contribution;
        linear_image[x, y] = color + (T{1} - coverage_fraction) * linear_image[x, y];
    });
}

template <std::floating_point T>
constexpr inline Vec2<T> texture_uv(const Vec3<T> &bary, const T v1_clipw, const T v2_clipw, const T v3_clipw,
                                    const Vec2<T> &uv1, const Vec2<T> &uv2, const Vec2<T> &uv3) noexcept
{
    T iwa = T{1} / v1_clipw;
    T iwb = T{1} / v2_clipw;
    T iwc = T{1} / v3_clipw;

    Vec2<T> uv1i = uv1 * iwa;
    Vec2<T> uv2i = uv2 * iwb;
    Vec2<T> uv3i = uv3 * iwc;

    Vec2<T> uv_interp = bary.x() * uv1i + bary.y() * uv2i + bary.z() * uv3i;
    T invw_interp = bary.x() * iwa + bary.y() * iwb + bary.z() * iwc;

    return uv_interp / invw_interp;
}

template <std::floating_point T>
constexpr inline Vec3<T> clip_to_screen_space(const Vec4<T> &clip_position, const Resolution &resolution) noexcept
{
    CHECK(abs(clip_position.w()) > T(1e-6));
    const Vec3<T> ndc = clip_position.xyz() / clip_position.w();
    const T sx = (ndc.x() * T{0.5} + T{0.5}) * static_cast<T>(resolution.width);
    const T sy = (T{1} - (ndc.y() * T{0.5} + T{0.5})) * static_cast<T>(resolution.height);
    return {sx, sy, ndc.z()};
}

template <std::floating_point T, VertexShader<T> VertexShader, FragmentShader<T> FragmentShader, size_t TileRows = 8,
          size_t TileCols = TileRows>
constexpr inline void render(const std::vector<Face> &faces, const std::vector<Vec3<T>> &positions,
                             const std::vector<Vec3<T>> &normals, const std::vector<Vec2<T>> &texture_coordinates,
                             const VertexShader &vertex_shader, const FragmentShader &fragment_shader, size_t aa_rows,
                             size_t aa_cols, IndexBuffer &index_buffer, std::vector<OutputVertex<T>> &vertex_buffer,
                             std::vector<Vec3<T>> &screen_coordinates, DepthBuffer<T> &depth_buffer,
                             ColorImage<T> &linear_image) noexcept
{
    using namespace std;
    CHECK(depth_buffer.width() == linear_image.width());
    CHECK(depth_buffer.height() == linear_image.height());

    const std::vector<Vec2<T>> sample_offsets = make_aa_grid<T>(aa_rows, aa_cols);

    const size_t vertex_count = Face::size * size(faces);
    vertex_buffer.resize(vertex_count);
    screen_coordinates.resize(vertex_count);
    if ((depth_buffer.width() != index_buffer.width()) || (depth_buffer.height() != index_buffer.height()) ||
        (depth_buffer.depth() != index_buffer.depth()))
    {
        index_buffer = IndexBuffer(depth_buffer.width(), depth_buffer.height(), depth_buffer.depth());
    }

    fill(execution::par_unseq, begin(index_buffer), end(index_buffer), numeric_limits<size_t>::max());
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(size(faces)), [&](size_t face_index) {
        const Face &face = faces[face_index];
        for (size_t vertex_index = 0; vertex_index < Face::size; ++vertex_index)
        {
            const size_t output_index = face_index * Face::size + vertex_index;
            const OutputVertex<T> output_vertex = vertex_shader({
                positions[face.vertex_indices[vertex_index]],
                normals[face.normal_indices[vertex_index]],
                texture_coordinates[face.texture_indices[vertex_index]],
            });

            vertex_buffer[output_index] = output_vertex;
            screen_coordinates[output_index] =
                clip_to_screen_space(output_vertex.clip_position, linear_image.resolution());
        }
    });
    rasterize_triangles<T, TileRows, TileCols>(screen_coordinates, sample_offsets, depth_buffer, index_buffer);
    render_vertices(screen_coordinates, vertex_buffer, fragment_shader, sample_offsets, index_buffer, linear_image);
}

template <std::floating_point T> class Renderer
{
  public:
    template <VertexShader<T> VertexShader, FragmentShader<T> FragmentShader, size_t TileRows = 8,
              size_t TileCols = TileRows>
    constexpr inline void render(const std::vector<Face> &faces, const std::vector<Vec3<T>> &positions,
                                 const std::vector<Vec3<T>> &normals, const std::vector<Vec2<T>> &texture_coordinates,
                                 const VertexShader &vertex_shader, const FragmentShader &fragment_shader,
                                 size_t aa_rows, size_t aa_cols, DepthBuffer<T> &depth_buffer,
                                 ColorImage<T> &linear_image) noexcept
    {
        ::render(faces, positions, normals, texture_coordinates, vertex_shader, fragment_shader, aa_rows, aa_cols,
                 index_buffer, vertex_buffer, screen_coordinates, depth_buffer, linear_image);
    }

  private:
    IndexBuffer index_buffer{};
    std::vector<OutputVertex<T>> vertex_buffer{};
    std::vector<Vec3<T>> screen_coordinates{};
};
