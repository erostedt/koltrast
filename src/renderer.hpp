#pragma once

#include <concepts>
#include <execution>
#include <vector>

#include "camera.hpp"
#include "counting_iterator.hpp"
#include "image.hpp"
#include "math.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
#include "shader.hpp"

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

template <std::floating_point T, FragmentShader<T> Shader, size_t AARows = 1, size_t AACols = AARows>
inline void render(ColorImage<T> &linear_image, const std::vector<OutputVertex<T>> &vertices,
                   const std::vector<Vec3<T>> &screen_coordinates, const Shader &fragment_shader,
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
                    const Vec3<T> sv1 = screen_coordinates[3 * index + 0];
                    const Vec3<T> sv2 = screen_coordinates[3 * index + 1];
                    const Vec3<T> sv3 = screen_coordinates[3 * index + 2];

                    const OutputVertex<T> v1 = vertices[3 * index + 0];
                    const OutputVertex<T> v2 = vertices[3 * index + 1];
                    const OutputVertex<T> v3 = vertices[3 * index + 2];

                    const Vec2<T> offset = offsets[ix, iy];

                    T px = (T)x + T(0.5) + offset.x();
                    T py = (T)y + T(0.5) + offset.y();

                    const Vec3<T> bary = barycentric({px, py}, sv1, sv2, sv3);
                    const Vec3<T> world_position =
                        bary.x() * v1.world_position + bary.y() * v2.world_position + bary.z() * v3.world_position;
                    const Vec3<T> world_normal =
                        bary.x() * v1.world_normal + bary.y() * v2.world_normal + bary.z() * v3.world_normal;
                    const Vec2<T> uv =
                        texture_uv(bary, v1.clip_position.w(), v2.clip_position.w(), v3.clip_position.w(),
                                   v1.texture_coordinates, v2.texture_coordinates, v3.texture_coordinates);

                    color = color + sample_contribution * fragment_shader({world_position, world_normal, uv});
                    ++coverage;
                }
            }
        }
        T coverage_fraction = (T)coverage * sample_contribution;
        linear_image[x, y] = color + (T{1} - coverage_fraction) * linear_image[x, y];
    });
}

// ADD AA templates
template <std::floating_point T> struct Renderer
{
    // TODO: Pass in generic vs/fs
    void render(const std::vector<Face> &faces, const std::vector<Vec3<T>> &positions,
                const std::vector<Vec3<T>> &normals, const std::vector<Vec2<T>> &texture_coordinates,
                const DefaultVertexShader<T> vs, const DefaultFragmentShader<T> &fs, DepthBuffer<T> depth_buffer,
                ColorImage<T> &linear_image)
    {
        CHECK(depth_buffer.resolution() == linear_image.resolution());

        using namespace std;
        const size_t vertex_count = Face::size * size(faces);
        vertex_buffer.resize(vertex_count);
        screen_coordinates.resize(vertex_count);
        const Resolution resolution = linear_image.resolution();
        if (resolution != index_buffer.resolution())
        {
            index_buffer = IndexBuffer<>(linear_image.width(), linear_image.height());
        }
        reset_index_buffer(index_buffer);
        for_each(execution::par_unseq, counting_iterator(0), counting_iterator(size(faces)), [&](size_t face_index) {
            const Face &face = faces[face_index];
            for (size_t vertex_index = 0; vertex_index < Face::size; ++vertex_index)
            {
                const size_t output_index = face_index * Face::size + vertex_index;
                const OutputVertex<T> output_vertex = vs({
                    positions[face.vertex_indices[vertex_index]],
                    normals[face.normal_indices[vertex_index]],
                    texture_coordinates[face.texture_indices[vertex_index]],
                });

                vertex_buffer[output_index] = output_vertex;
                screen_coordinates[output_index] = clip_to_screen_space(output_vertex.clip_position, resolution);
            }
        });
        rasterize_triangles(screen_coordinates, depth_buffer, index_buffer);
        ::render(linear_image, vertex_buffer, screen_coordinates, fs, index_buffer);
    }

    IndexBuffer<> index_buffer;
    std::vector<OutputVertex<T>> vertex_buffer;
    std::vector<Vec3<T>> screen_coordinates;
};

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
