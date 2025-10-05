#pragma once

#include <algorithm>
#include <concepts>

#include "color.hpp"
#include "types.hpp"

template <std::floating_point T> struct PointLight
{
    Vec3<T> position;
    RGBA<T> color;
    T specular;
};

template <std::floating_point T> struct DirectionalLight
{
    Vec3<T> direction;
    RGBA<T> color;
    T specular;
};

template <std::floating_point T>
[[nodiscard]] constexpr inline RGBA<T> directional_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                         const Vec3<T> &camera_position, T object_shininess,
                                                         const Vec3<T> &light_direction, const RGBA<T> &light_color,
                                                         T specular) noexcept
{
    const Vec3<T> norm = *world_normal.normalized();
    const Vec3<T> view_direction = *(camera_position - world_position).normalized();

    // Diffuse
    const T diffuse = std::max(norm.dot(light_direction), T{0});
    const RGBA<T> diffuse_color = diffuse * light_color;

    // Specular
    const Vec3<T> reflection_direction = *(norm * (T{2} * norm.dot(light_direction)) - light_direction).normalized();
    const T spec = specular * std::pow(std::max(view_direction.dot(reflection_direction), T{0}), object_shininess);
    const RGBA<T> specular_color = spec * light_color;

    return diffuse_color + specular_color;
}

template <std::floating_point T>
[[nodiscard]] constexpr inline RGBA<T> sample_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                    const Vec3<T> &camera_position, T object_shininess,
                                                    const PointLight<T> &pl) noexcept
{
    const Vec3<T> light_direction = *(pl.position - world_position).normalized();
    return directional_light(world_position, world_normal, camera_position, object_shininess, light_direction, pl.color,
                             pl.specular);
}

template <std::floating_point T>
[[nodiscard]] constexpr inline RGBA<T> sample_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                    const Vec3<T> &camera_position, T object_shininess,
                                                    const DirectionalLight<T> &dl) noexcept
{
    return directional_light(world_position, world_normal, camera_position, object_shininess, dl.direction, dl.color,
                             dl.specular);
}
