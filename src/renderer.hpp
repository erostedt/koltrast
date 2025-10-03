#pragma once

#include "counting_iterator.hpp"
#include "light.hpp"
#include "math.hpp"
#include "rasterizer.hpp"
#include "texture.hpp"
#include <concepts>

template <std::floating_point T> struct Fragment
{
    Vec3<T> world_position;
    Vec3<T> world_normal;
    Vec2<T> uv;
};

template <typename F, typename T>
concept FragmentShader = std::floating_point<T> && requires(F f, const Fragment<T> &fragment) {
    { f(fragment) } -> std::same_as<RGB<T>>;
};

template <std::floating_point T> constexpr inline RGB<T> srgb_to_linear(const RGB<T> &c) noexcept
{
    auto convert = [](T v) {
        if (v <= T(0.04045))
            return v / T(12.92);
        else
            return std::pow((v + T(0.055)) / T(1.055), T(2.4));
    };
    return {convert(c.r), convert(c.g), convert(c.b)};
}

template <std::floating_point T> constexpr inline RGB<T> linear_to_srgb(const RGB<T> &c) noexcept
{
    auto convert = [](T v) {
        if (v <= T(0.0031308))
            return v * T(12.92);
        else
            return T(1.055) * std::pow(v, T(1.0 / 2.4)) - T(0.055);
    };
    return {convert(c.r), convert(c.g), convert(c.b)};
}

template <std::floating_point T>
constexpr inline RGB<T> sample_bilinear(const Vec2<T> &uv, const Texture<T> &texture) noexcept
{
    T x = std::clamp(uv.x() * (T)(texture.width() - 1), T{0}, (T)(texture.width() - 1));
    T y = std::clamp(uv.y() * (T)(texture.height() - 1), T{0}, (T)(texture.height() - 1));

    size_t fx = floor_to_size(x);
    size_t cx = ceil_to_size(x);
    size_t fy = floor_to_size(y);
    size_t cy = ceil_to_size(y);

    RGB<T> a = srgb_to_linear(texture[fx, fy]);
    RGB<T> b = srgb_to_linear(texture[cx, fy]);
    RGB<T> c = srgb_to_linear(texture[fx, cy]);
    RGB<T> d = srgb_to_linear(texture[cx, cy]);

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
    return srgb_to_linear(texture[tx, ty]);
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

template <std::floating_point T> struct DefaultShader
{
    Vec3<T> camera_position;
    Lights<T> lights;
    Texture<T> texture;
    T object_shininess;

    [[nodiscard]] constexpr inline RGB<T> operator()(const Fragment<T> &fragment) const
    {
        const RGB<T> object_color = sample_bilinear(fragment.uv, texture);
        const RGB<T> light =
            accumulate_light(lights, fragment.world_position, fragment.world_normal, camera_position, object_shininess);
        return light * object_color;
    }
};

template <std::floating_point T, FragmentShader<T> Shader, size_t AARows = 1, size_t AACols = AARows>
inline void render(ColorImage<T> &linear_image, const std::vector<Face> &faces, const VertexBuffer<T> &vertex_buffer,
                   const std::vector<Vec2<T>> &texture_coordinates, const Shader &fragment_shader,
                   const IndexBuffer<AARows, AACols> &index_buffer) noexcept
{
    using namespace std;
    static constexpr Matrix<Vec2<T>, AARows, AACols> offsets = make_aa_grid<T, AARows, AACols>();

    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();
        const Matrix<size_t, AARows, AACols> &indicies = index_buffer[x, y];
        if (all_of(begin(indicies), end(indicies), [](size_t idx) { return idx == numeric_limits<size_t>::max(); }))
        {
            return;
        }

        const T sample_contribution = T{1} / (AARows * AACols);

        size_t coverage = 0;
        RGB<T> color = BLACK<T>;
        for (size_t iy = 0; iy < AARows; ++iy)
        {
            for (size_t ix = 0; ix < AACols; ++ix)
            {
                size_t index = indicies[ix, iy];
                if (index != numeric_limits<size_t>::max())
                {
                    const Face &face = faces[index];
                    const Vec4<T> sv1 = vertex_buffer.screen_coordinates[face.vertex_indices[0]];
                    const Vec4<T> sv2 = vertex_buffer.screen_coordinates[face.vertex_indices[1]];
                    const Vec4<T> sv3 = vertex_buffer.screen_coordinates[face.vertex_indices[2]];

                    const Vec3<T> wv1 = vertex_buffer.positions[face.vertex_indices[0]];
                    const Vec3<T> wv2 = vertex_buffer.positions[face.vertex_indices[1]];
                    const Vec3<T> wv3 = vertex_buffer.positions[face.vertex_indices[2]];

                    const Vec3<T> wn1 = vertex_buffer.normals[face.normal_indices[0]];
                    const Vec3<T> wn2 = vertex_buffer.normals[face.normal_indices[1]];
                    const Vec3<T> wn3 = vertex_buffer.normals[face.normal_indices[2]];

                    const Vec2<T> uv1 = texture_coordinates[face.texture_indices[0]];
                    const Vec2<T> uv2 = texture_coordinates[face.texture_indices[1]];
                    const Vec2<T> uv3 = texture_coordinates[face.texture_indices[2]];

                    const Vec2<T> offset = offsets[ix, iy];

                    T px = (T)x + T(0.5) + offset.x();
                    T py = (T)y + T(0.5) + offset.y();

                    const Vec3<T> bary = barycentric({px, py}, sv1, sv2, sv3);
                    const Vec3<T> world_position = bary.x() * wv1 + bary.y() * wv2 + bary.z() * wv3;
                    const Vec3<T> world_normal = bary.x() * wn1 + bary.y() * wn2 + bary.z() * wn3;
                    const Vec2<T> uv = texture_uv(bary, sv1, sv2, sv3, uv1, uv2, uv3);

                    color = color + sample_contribution * fragment_shader({world_position, world_normal, uv});
                    ++coverage;
                }
            }
        }
        T coverage_fraction = (T)coverage * sample_contribution;
        linear_image[x, y] = color + (T{1} - coverage_fraction) * linear_image[x, y];
    });
}

template <std::floating_point T> inline void linear_to_srgb(ColorImage<T> &image) noexcept
{
    using namespace std;
    std::transform(execution::par, begin(image), end(image), begin(image),
                   [](const RGB<T> &color) { return linear_to_srgb<T>(color); });
}

template <std::floating_point T> RGB<T> sample_cubemap(const Vec3<T> &direction, const Image<RGB<T>> &cubemap)
{
    T phi = std::atan2(direction.z(), direction.x());
    T theta = std::acos(direction.y());

    T pi = std::numbers::pi_v<f32>;
    T u = (phi + pi) / (T{2} * pi);
    T v = theta / pi;

    return sample_bilinear({u, v}, cubemap);
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
