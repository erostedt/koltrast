#pragma once

#include <concepts>
#include <numeric>
#include <vector>

#include "image.hpp"
#include "light.hpp"
#include "math.hpp"
#include "texture.hpp"
#include "types.hpp"

template <std::floating_point T> struct InputVertex
{
    Vec3<T> position;
    Vec3<T> normal;
    Vec2<T> texture_coordinates;
};

template <std::floating_point T> struct OutputVertex
{
    Vec3<T> position;
    Vec3<T> normal;
    Vec2<T> texture_coordinates;
    Vec4<T> clip_position;
};

template <typename F, typename T>
concept VertexShader = std::floating_point<T> && requires(F f, const InputVertex<T> &vertex) {
    { f(vertex) } -> std::same_as<OutputVertex<T>>;
};

template <std::floating_point T> class DefaultVertexShader
{
  public:
    DefaultVertexShader(const Mat4x4<T> &model_matrix, const Mat4x4<T> &view_matrix, const Mat4x4<T> &projection_matrix)
        : _model_matrix(model_matrix), _view_projection_matrix(projection_matrix * view_matrix)
    {
        Mat3x3<T> mat3 = {model_matrix[0, 0], model_matrix[0, 1], model_matrix[0, 2],
                          model_matrix[1, 0], model_matrix[1, 1], model_matrix[1, 2],
                          model_matrix[2, 0], model_matrix[2, 1], model_matrix[2, 2]};

        _normal_matrix = mat3.inverse().value().transposed();
    }

    [[nodiscard]] constexpr inline OutputVertex<T> operator()(const InputVertex<T> &vertex) const
    {
        using namespace std;
        const Vec4<T> world_position =
            _model_matrix * Vec4<T>{vertex.position.x(), vertex.position.y(), vertex.position.z(), T{1}};
        const T invw = T{1} / world_position.w();

        OutputVertex<T> out;
        out.position.x() = world_position.x() * invw;
        out.position.y() = world_position.y() * invw;
        out.position.z() = world_position.z() * invw;
        out.normal = *(_normal_matrix * vertex.normal).normalized();
        out.texture_coordinates = vertex.texture_coordinates;
        out.clip_position = _view_projection_matrix * world_position;
        return out;
    }

  private:
    Mat4x4<T> _model_matrix;
    Mat4x4<T> _view_projection_matrix;
    Mat3x3<T> _normal_matrix;
};

// TODO: (eric) merge with InputVertex?
template <std::floating_point T> struct Fragment
{
    Vec3<T> position;
    Vec3<T> normal;
    Vec2<T> uv;
};

template <typename F, typename T>
concept FragmentShader = std::floating_point<T> && requires(F f, const Fragment<T> &fragment) {
    { f(fragment) } -> std::same_as<RGBA<T>>;
};

template <typename F, typename T, typename Shader>
concept AAFunction = std::floating_point<T> && FragmentShader<Shader, T> &&
                     requires(const F &f, const Shader &fragment_shader, const std::vector<Fragment<T>> &fragments,
                              size_t sample_count) {
                         { f(fragment_shader, fragments, sample_count) } -> std::same_as<RGBA<T>>;
                     };

template <std::floating_point T, FragmentShader<T> FragmentShader> struct NoAA
{
    [[nodiscard]] constexpr inline RGBA<T> operator()(const FragmentShader &fragment_shader,
                                                      const std::vector<Fragment<T>> &fragments,
                                                      size_t /*sample_count*/) const noexcept
    {
        return fragment_shader(fragments.front());
    }
};

template <std::floating_point T, FragmentShader<T> FragmentShader> struct SSAA
{
    [[nodiscard]] constexpr inline RGBA<T> operator()(const FragmentShader &fragment_shader,
                                                      const std::vector<Fragment<T>> &fragments,
                                                      size_t sample_count) const noexcept
    {
        const T sample_contribution = T{1} / T(sample_count);
        using namespace std;
        return transform_reduce(begin(fragments), end(fragments), BLACK<T>, std::plus<RGBA<T>>{},
                                [&](const Fragment<T> &fragment) {
                                    const RGBA<T> color = fragment_shader(fragment);
                                    return sample_contribution * color;
                                });
    }
};

template <std::floating_point T, FragmentShader<T> FragmentShader> struct MSAA
{
    [[nodiscard]] constexpr inline RGBA<T> operator()(const FragmentShader &fragment_shader,
                                                      const std::vector<Fragment<T>> &fragments,
                                                      size_t sample_count) const noexcept
    {
        using namespace std;
        CHECK(size(fragments) > 0);
        const T sample_contibution = T{1} / T(size(fragments));
        const T coverage = T(size(fragments)) / T(sample_count);
        Vec3<T> position{};
        Vec3<T> normal{};
        Vec2<T> uv{};
        for (const Fragment<T> &fragment : fragments)
        {
            position += fragment.position;
            normal += fragment.normal;
            uv += fragment.uv;
        }
        position *= sample_contibution;
        normal *= sample_contibution;
        uv *= sample_contibution;
        return coverage * fragment_shader({position, normal, uv});
    }
};

template <std::floating_point T, typename ColorType>
constexpr inline RGBA<T> sample_nearest_neighbor(const Vec2<T> &uv, const Image<RGBA<ColorType>> &texture) noexcept
{
    const size_t tx = floor_to_size(uv.x() * (T)(texture.width() - 1));
    const size_t ty = floor_to_size(uv.y() * (T)(texture.height() - 1));
    return srgb_to_linear<T>(texture[tx, ty]);
}

template <std::floating_point T, typename ColorType>
constexpr inline RGBA<T> sample_bilinear(const Vec2<T> &uv, const Image<RGBA<ColorType>> &texture) noexcept
{
    const T x = std::clamp(uv.x() * (T)(texture.width() - 1), T{0}, (T)(texture.width() - 1));
    const T y = std::clamp(uv.y() * (T)(texture.height() - 1), T{0}, (T)(texture.height() - 1));

    const size_t fx = floor_to_size(x);
    const size_t cx = ceil_to_size(x);
    const size_t fy = floor_to_size(y);
    const size_t cy = ceil_to_size(y);

    const RGBA<T> a = srgb_to_linear<T>(texture[fx, fy]);
    const RGBA<T> b = srgb_to_linear<T>(texture[cx, fy]);
    const RGBA<T> c = srgb_to_linear<T>(texture[fx, cy]);
    const RGBA<T> d = srgb_to_linear<T>(texture[cx, cy]);

    const T s = x - (T)fx;
    const T t = y - (T)fy;

    return (T{1} - s) * (T{1} - t) * a + s * (T{1} - t) * b + (T{1} - s) * t * c + s * t * d;
}

template <std::floating_point T> RGBA<T> sample_cubemap(const Vec3<T> &direction, const Image<RGBA<T>> &cubemap)
{
    T phi = std::atan2(direction.z(), direction.x());
    T theta = std::acos(direction.y());

    T pi = std::numbers::pi_v<f32>;
    T u = (phi + pi) / (T{2} * pi);
    T v = theta / pi;

    return sample_bilinear<T>({u, v}, cubemap);
}

template <std::floating_point T> struct DefaultFragmentShader
{
    Vec3<T> camera_position;
    Texture texture;
    T object_shininess;
    std::vector<PointLight<T>> point_lights;
    std::vector<DirectionalLight<T>> directional_lights;
    T ambient = T{0};

    [[nodiscard]] constexpr inline RGBA<T> operator()(const Fragment<T> &fragment) const
    {
        const RGBA<T> object_color = sample_bilinear(fragment.uv, texture);

        RGB<T> total = {ambient, ambient, ambient};
        for (const auto &pl : point_lights)
        {
            const RGB<T> light =
                sample_light(fragment.position, fragment.normal, camera_position, object_shininess, pl);
            total += light;
        }

        for (const auto &dl : directional_lights)
        {
            const RGB<T> light =
                sample_light(fragment.position, fragment.normal, camera_position, object_shininess, dl);
            total += light;
        }
        total = std::clamp(total, T{0}, T{1});
        const RGBA<T> light = total.with_alpha(T{1});
        return light * object_color;
    }
};
