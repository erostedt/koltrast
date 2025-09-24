#pragma once

#include "counting_iterator.hpp"
#include "rasterizer.hpp"
#include "texture.hpp"

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

constexpr inline Vec3f barycentric(const Vec2f &at, const Vec4f &v1, const Vec4f &v2, const Vec4f &v3) noexcept
{
    const Vec4f p = {at.x(), at.y(), 0.0f, 0.0f};
    f32 w = 1.0f / edge_function(v1, v2, v3);
    f32 a = edge_function(v2, v3, p) * w;
    f32 b = edge_function(v3, v1, p) * w;
    return {a, b, 1.0f - a - b};
}

constexpr inline Vec2f interpolate_uv(const Vec3f &bary, const Vec4f &v1, const Vec4f &v2, const Vec4f &v3,
                                      const Vec2f &uv1, const Vec2f &uv2, const Vec2f &uv3) noexcept
{
    f32 wa = v1.w();
    f32 wb = v2.w();
    f32 wc = v3.w();

    f32 iwa = 1.0f / wa;
    f32 iwb = 1.0f / wb;
    f32 iwc = 1.0f / wc;

    Vec2f uv1i = uv1 * iwa;
    Vec2f uv2i = uv2 * iwb;
    Vec2f uv3i = uv3 * iwc;

    Vec2f uv_interp = bary.x() * uv1i + bary.y() * uv2i + bary.z() * uv3i;
    f32 invw_interp = bary.x() * iwa + bary.y() * iwb + bary.z() * iwc;

    return uv_interp / invw_interp;
}

inline void render_triangles(ColorImage &image, const std::vector<RGB<f32>> &colors,
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

inline void render_triangles(ColorImage &image, const std::vector<Face> &faces,
                             const std::vector<Vec4f> &screen_vertices, const std::vector<Vec2f> &texture_coordinates,
                             const Texture &texture, const IndexBuffer &index_buffer) noexcept
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

            const auto bary = barycentric({(f32)x + 0.5f, (f32)y + 0.5f}, v1, v2, v3);
            const auto uv = interpolate_uv(bary, v1, v2, v3, uv1, uv2, uv3);
            image[x, y] = sample_bilinear(uv, texture);
        }
    });
}
