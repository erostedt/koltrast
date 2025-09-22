#include "camera.hpp"

using Vec3f = Vector<f32, 3>;
using Vec2f = Vector<f32, 2>;

inline RGB<f32> point_light(const Vec3f &world_position, const Vec3f &world_normal, const Vec3f &light_position,
                            const Vec3f &camera_position, const RGB<f32> &object_color, const RGB<f32> &light_color,
                            float ambient, float specular, float shininess)
{
    const Vec3f norm = *world_normal.normalized();
    const Vec3f light_direction = *(light_position - world_position).normalized();
    const Vec3f view_direction = *(camera_position - world_position).normalized();

    // Ambient
    const RGB<f32> ambient_color = {ambient * light_color.r, ambient * light_color.g, ambient * light_color.b};

    // Diffuse
    const float diffuse = std::max(norm.dot(light_direction), 0.0f);
    const RGB<f32> diffuse_color = {diffuse * light_color.r, diffuse * light_color.g, diffuse * light_color.b};

    // Specular
    const Vec3f reflectDir = *(norm * (2.0f * norm.dot(light_direction)) - light_direction).normalized();
    const float spec = std::pow(std::max(view_direction.dot(reflectDir), 0.0f), shininess);
    const RGB<f32> specular_color = {specular * spec * light_color.r, specular * spec * light_color.g,
                                     specular * spec * light_color.b};

    return {
        (ambient_color.r + diffuse_color.r + specular_color.r) * object_color.r,
        (ambient_color.g + diffuse_color.g + specular_color.g) * object_color.g,
        (ambient_color.b + diffuse_color.b + specular_color.b) * object_color.b,
    };
}
