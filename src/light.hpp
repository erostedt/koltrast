#include "camera.hpp"

using Vec3f = Vector<f32, 3>;
using Vec2f = Vector<f32, 2>;

struct PointLight
{
    Vec3f position;
    RGB<f32> color;
    f32 specular;
};

struct DirectionalLight
{
    Vec3f direction;
    RGB<f32> color;
    f32 specular;
};

[[nodiscard]] constexpr inline RGB<f32> directional_light(const Vec3f &world_position, const Vec3f &world_normal,
                                                          const Vec3f &camera_position, f32 object_shininess,
                                                          const Vec3f &light_direction, const RGB<f32> &light_color,
                                                          f32 specular) noexcept
{
    const Vec3f norm = *world_normal.normalized();
    const Vec3f view_direction = *(camera_position - world_position).normalized();

    // Diffuse
    const float diffuse = std::max(norm.dot(light_direction), 0.0f);
    const RGB<f32> diffuse_color = {diffuse * light_color.r, diffuse * light_color.g, diffuse * light_color.b};

    // Specular
    const Vec3f reflection_direction = *(norm * (2.0f * norm.dot(light_direction)) - light_direction).normalized();
    const float spec = std::pow(std::max(view_direction.dot(reflection_direction), 0.0f), object_shininess);
    const RGB<f32> specular_color = {specular * spec * light_color.r, specular * spec * light_color.g,
                                     specular * spec * light_color.b};

    return {
        (diffuse_color.r + specular_color.r),
        (diffuse_color.g + specular_color.g),
        (diffuse_color.b + specular_color.b),
    };
}

[[nodiscard]] constexpr inline RGB<f32> sample_light(const Vec3f &world_position, const Vec3f &world_normal,
                                                     const Vec3f &camera_position, f32 object_shininess,
                                                     const PointLight &pl) noexcept
{
    const Vec3f light_direction = *(pl.position - world_position).normalized();
    return directional_light(world_position, world_normal, camera_position, object_shininess, light_direction, pl.color,
                             pl.specular);
}

[[nodiscard]] constexpr inline RGB<f32> sample_light(const Vec3f &world_position, const Vec3f &world_normal,
                                                     const Vec3f &camera_position, f32 object_shininess,
                                                     const DirectionalLight &dl) noexcept
{
    return directional_light(world_position, world_normal, camera_position, object_shininess, dl.direction, dl.color,
                             dl.specular);
}
