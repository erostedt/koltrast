#pragma once

#include "counting_iterator.hpp"
#include "light.hpp"
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
constexpr inline Vec2<T> texture_uv(const Vec3<T> &bary, const Vec4<T> &v1, const Vec4<T> &v2, const Vec4<T> &v3,
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
            const auto uv = texture_uv(bary, v1, v2, v3, uv1, uv2, uv3);
            image[x, y] = sample_bilinear(uv, texture);
        }
    });
}

template <std::floating_point T>
void apply_lighting(ColorImage<T> &image, const Vec3<T> &camera_position, const std::vector<Face> &faces,
                    const std::vector<Vec4<T>> &screen_vertices, const std::vector<Vec4<T>> &world_vertices,
                    const std::vector<Vec3<T>> &world_normals, const IndexBuffer &index_buffer)
{
    PointLight<T> point_light{
        .position = {T(0), T(1), T(2)},
        .color = {T{1}, T{1}, T{1}},
        .specular = T(0.8),
    };

    T ambient = T(0.3);
    T shininess = T(16);
    using namespace std;
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();
        size_t index = index_buffer[x, y];
        if (index != std::numeric_limits<size_t>::max())
        {
            const auto &face = faces[index];
            const auto sv1 = screen_vertices[face.vertex_indices[0]];
            const auto sv2 = screen_vertices[face.vertex_indices[1]];
            const auto sv3 = screen_vertices[face.vertex_indices[2]];

            const auto wv1 = world_vertices[face.vertex_indices[0]];
            const auto wv2 = world_vertices[face.vertex_indices[1]];
            const auto wv3 = world_vertices[face.vertex_indices[2]];

            const auto wn1 = world_normals[face.normal_indices[0]];
            const auto wn2 = world_normals[face.normal_indices[1]];
            const auto wn3 = world_normals[face.normal_indices[2]];

            const auto bary = barycentric({(T)x + T(0.5), (T)y + T(0.5)}, sv1, sv2, sv3);
            const auto world_position = (bary.x() * wv1 + bary.y() * wv2 + bary.z() * wv3).xyz();
            const auto world_normal = (bary.x() * wn1 + bary.y() * wn2 + bary.z() * wn3);

            const RGB<T> light = sample_light(world_position, world_normal, camera_position, shininess, point_light);

            const RGB<T> object_color = image[x, y];
            image[x, y] = {(ambient + light.r) * object_color.r, (ambient + light.g) * object_color.g,
                           (ambient + light.b) * object_color.b};
        }
    });
}

template <std::floating_point T> RGB<T> sample_cubemap(const Vec3<T> &direction, const Image<RGB<T>> &cubemap)
{
    T phi = std::atan2(direction.z(), direction.x());
    T theta = std::acos(direction.y());

    T pi = std::numbers::pi_v<f32>;
    T u = (phi + pi) / (T{2} * pi);
    T v = theta / pi;

    // TODO: (ecrt) Use bilinear interpolation
    size_t px = (size_t)(u * (T)(cubemap.width() - 1));
    size_t py = (size_t)(v * (T)(cubemap.height() - 1));

    return cubemap[px, py];
}

template <std::floating_point T> struct ViewPort
{
    Vec3<T> pixel_delta_u;
    Vec3<T> pixel_delta_v;
    Vec3<T> pixel_top_left;
};

template <std::floating_point T>
ViewPort<T> view_port(const Camera<T> &camera, const Vec3<T> &position, const Mat4x4<T> &view)
{
    T aspect_ratio = (T)camera.resolution.width / (T)camera.resolution.height;

    T angle = radians(camera.vertical_field_of_view);
    T height = std::tan(angle / T{2});
    T viewport_height = T{2} * height * camera.focal_length;
    T viewport_width = viewport_height * aspect_ratio;

    Vec3<T> right = {view[0, 0], view[0, 1], view[0, 2]};
    Vec3<T> up = {view[1, 0], view[1, 1], view[1, 2]};
    Vec3<T> forward = {view[2, 0], view[2, 1], view[2, 2]};

    Vec3<T> viewport_u = viewport_width * right;
    Vec3<T> viewport_v = -viewport_height * up;

    Vec3<T> viewport_center = position - camera.focal_length * forward;
    Vec3<T> viewport_upper_left = viewport_center - viewport_u / T{2} - viewport_v / T{2};

    Vec3<T> pixel_delta_u = viewport_u / (T)camera.resolution.width;
    Vec3<T> pixel_delta_v = viewport_v / (T)camera.resolution.height;

    Vec3<T> pixel_top_left = viewport_upper_left + T(0.5) * (pixel_delta_u + pixel_delta_v);
    return {pixel_delta_u, pixel_delta_v, pixel_top_left};
}
