#pragma once
#include "camera.hpp"
#include "types.hpp"
#include <algorithm>
#include <concepts>
#include <vector>

template <std::floating_point T> struct PointLight
{
    Vec3<T> position;
    RGB<T> color;
    T specular;
};

template <std::floating_point T> struct DirectionalLight
{
    Vec3<T> direction;
    RGB<T> color;
    T specular;
};

template <std::floating_point T> struct Lights
{
    T ambient = T{0};
    std::vector<PointLight<T>> point_lights;
    std::vector<DirectionalLight<T>> directional_lights;
};

template <std::floating_point T>
RGB<T> accumulate_light(const Lights<T> &lights, const Vec3<T> &world_position, const Vec3<T> &world_normal,
                        const Vec3<T> &camera_position, T object_shininess)
{
    RGB<T> total = {lights.ambient, lights.ambient, lights.ambient};
    for (const auto &pl : lights.point_lights)
    {
        const RGB<T> light = sample_light(world_position, world_normal, camera_position, object_shininess, pl);
        total.r += light.r;
        total.g += light.g;
        total.b += light.b;
    }

    for (const auto &dl : lights.directional_lights)
    {
        const RGB<T> light = sample_light(world_position, world_normal, camera_position, object_shininess, dl);
        total.r += light.r;
        total.g += light.g;
        total.b += light.b;
    }
    T r = std::clamp(total.r, T{0}, T{1});
    T g = std::clamp(total.g, T{0}, T{1});
    T b = std::clamp(total.b, T{0}, T{1});
    return {r, g, b};
}

template <std::floating_point T>
[[nodiscard]] constexpr inline RGB<T> directional_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                        const Vec3<T> &camera_position, T object_shininess,
                                                        const Vec3<T> &light_direction, const RGB<T> &light_color,
                                                        T specular) noexcept
{
    const Vec3<T> norm = *world_normal.normalized();
    const Vec3<T> view_direction = *(camera_position - world_position).normalized();

    // Diffuse
    const T diffuse = std::max(norm.dot(light_direction), T{0});
    const RGB<T> diffuse_color = {diffuse * light_color.r, diffuse * light_color.g, diffuse * light_color.b};

    // Specular
    const Vec3<T> reflection_direction = *(norm * (T{2} * norm.dot(light_direction)) - light_direction).normalized();
    const T spec = std::pow(std::max(view_direction.dot(reflection_direction), T{0}), object_shininess);
    const RGB<T> specular_color = {specular * spec * light_color.r, specular * spec * light_color.g,
                                   specular * spec * light_color.b};

    return {
        (diffuse_color.r + specular_color.r),
        (diffuse_color.g + specular_color.g),
        (diffuse_color.b + specular_color.b),
    };
}

template <std::floating_point T>
[[nodiscard]] constexpr inline RGB<T> sample_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                   const Vec3<T> &camera_position, T object_shininess,
                                                   const PointLight<T> &pl) noexcept
{
    const Vec3<T> light_direction = *(pl.position - world_position).normalized();
    return directional_light(world_position, world_normal, camera_position, object_shininess, light_direction, pl.color,
                             pl.specular);
}

template <std::floating_point T>
[[nodiscard]] constexpr inline RGB<T> sample_light(const Vec3<T> &world_position, const Vec3<T> &world_normal,
                                                   const Vec3<T> &camera_position, T object_shininess,
                                                   const DirectionalLight<T> &dl) noexcept
{
    return directional_light(world_position, world_normal, camera_position, object_shininess, dl.direction, dl.color,
                             dl.specular);
}
