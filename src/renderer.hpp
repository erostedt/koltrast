#pragma once

#include "counting_iterator.hpp"
#include "rasterizer.hpp"
#include "texture.hpp"
#include <concepts>

template <std::floating_point T>
constexpr inline RGB<T> sample_bilinear(const Vec2<T> &uv, const Texture<T> &texture) noexcept
{
    T x = std::clamp(uv.x() * (T)(texture.width() - 1), T{0}, (T)(texture.width() - 1));
    T y = std::clamp(uv.y() * (T)(texture.height() - 1), T{0}, (T)(texture.height() - 1));

    size_t fx = floor_to_size(x);
    size_t cx = ceil_to_size(x);
    size_t fy = floor_to_size(y);
    size_t cy = ceil_to_size(y);

    RGB<T> a = texture[fx, fy];
    RGB<T> b = texture[cx, fy];
    RGB<T> c = texture[fx, cy];
    RGB<T> d = texture[cx, cy];

    T s = x - (T)fx;
    T t = y - (T)fy;

    T red = (T{1} - s) * (T{1} - t) * a.r + s * (T{1} - t) * b.r + (T{1} - s) * t * c.r + s * t * d.r;
    T green = (T{1} - s) * (T{1} - t) * a.g + s * (T{1} - t) * b.g + (T{1} - s) * t * c.g + s * t * d.g;
    T blue = (T{1} - s) * (T{1} - t) * a.b + s * (T{1} - t) * b.b + (T{1} - s) * t * c.b + s * t * d.b;
    return {red, green, blue};
}

template <std::floating_point T>
constexpr inline RGB<T> sample_nearest_neighbor(const Vec2<T> &uv, const Texture<T> &texture) noexcept
{
    size_t tx = floor_to_size(uv.x() * (T)(texture.width() - 1));
    size_t ty = floor_to_size(uv.y() * (T)(texture.height() - 1));
    return texture[tx, ty];
}

template <std::floating_point T>
constexpr inline Vec3<T> barycentric(const Vec2<T> &at, const Vec4<T> &v1, const Vec4<T> &v2,
                                     const Vec4<T> &v3) noexcept
{
    const Vec4<T> p = {at.x(), at.y(), T{0}, T{0}};
    T w = T{1} / edge_function(v1, v2, v3);
    T a = edge_function(v2, v3, p) * w;
    T b = edge_function(v3, v1, p) * w;
    return {a, b, T{1} - a - b};
}

template <std::floating_point T>
constexpr inline Vec2<T> interpolate_uv(const Vec3<T> &bary, const Vec4<T> &v1, const Vec4<T> &v2, const Vec4<T> &v3,
                                        const Vec2<T> &uv1, const Vec2<T> &uv2, const Vec2<T> &uv3) noexcept
{
    T wa = v1.w();
    T wb = v2.w();
    T wc = v3.w();

    T iwa = T{1} / wa;
    T iwb = T{1} / wb;
    T iwc = T{1} / wc;

    Vec2<T> uv1i = uv1 * iwa;
    Vec2<T> uv2i = uv2 * iwb;
    Vec2<T> uv3i = uv3 * iwc;

    Vec2<T> uv_interp = bary.x() * uv1i + bary.y() * uv2i + bary.z() * uv3i;
    T invw_interp = bary.x() * iwa + bary.y() * iwb + bary.z() * iwc;

    return uv_interp / invw_interp;
}

template <std::floating_point T>
inline void render_triangles(ColorImage<T> &image, const std::vector<RGB<T>> &colors,
                             const IndexBuffer &index_buffer) noexcept
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

template <std::floating_point T>
inline void render_triangles(ColorImage<T> &image, const std::vector<Face> &faces,
                             const std::vector<Vec4<T>> &screen_vertices,
                             const std::vector<Vec2<T>> &texture_coordinates, const Texture<T> &texture,
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

            const auto bary = barycentric({(T)x + T(0.5), (T)y + T(0.5)}, v1, v2, v3);
            const auto uv = interpolate_uv(bary, v1, v2, v3, uv1, uv2, uv3);
            image[x, y] = sample_bilinear(uv, texture);
        }
    });
}
