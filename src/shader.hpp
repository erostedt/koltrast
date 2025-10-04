#pragma once

#include <concepts>

#include "image.hpp"
#include "light.hpp"
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
    // TODO: Remove world prefix
    Vec3<T> world_position;
    Vec3<T> world_normal;
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
        out.world_position.x() = world_position.x() * invw;
        out.world_position.y() = world_position.y() * invw;
        out.world_position.z() = world_position.z() * invw;
        out.world_normal = *(_normal_matrix * vertex.normal).normalized();
        out.texture_coordinates = vertex.texture_coordinates;
        out.clip_position = _view_projection_matrix * world_position;
        return out;
    }

  private:
    Mat4x4<T> _model_matrix;
    Mat4x4<T> _view_projection_matrix;
    Mat3x3<T> _normal_matrix;
};

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

template <std::floating_point T> struct DefaultFragmentShader
{
    Vec3<T> camera_position;
    Lights<T> lights;
    Texture texture;
    T object_shininess;

    [[nodiscard]] constexpr inline RGB<T> operator()(const Fragment<T> &fragment) const
    {
        const RGB<T> object_color = sample_bilinear(fragment.uv, texture);
        const RGB<T> light =
            accumulate_light(lights, fragment.world_position, fragment.world_normal, camera_position, object_shininess);
        return light * object_color;
    }
};
