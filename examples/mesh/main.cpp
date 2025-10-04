#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "camera.hpp"
#include "image.hpp"
#include "obj.hpp"
#include "renderer.hpp"
#include "texture.hpp"
#include "transform.hpp"
#include "types.hpp"

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << "./mesh obj_file texture_file > output.ppm";
        return 1;
    }

    const auto mesh = load_obj<f32>(argv[1]);
    const auto texture = load_texture(argv[2]);

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const Vec3f camera_position = {0.0f, 0.0f, 2.0f};
    const auto model = model_matrix<f32>({0.0f, 0.0f, 0.0f}, {0.0f, 45.0f, 0.0f}, {1.0f, 1.0f, 1.0f});
    const auto view = look_at<f32>(camera_position, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    ColorImage<f32> image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer<f32>(camera.resolution.width, camera.resolution.height);

    Renderer<f32> renderer;

    const DefaultFragmentShader<f32> fragment_shader = {
        .camera_position = camera_position,
        .texture = texture,
        .object_shininess = 16.0f,
        .point_lights = {{.position = {0.0f, 1.0f, 2.0f}, .color = {1.0f, 1.0f, 1.0f}, .specular = 0.8f}},
        .directional_lights = {},
        .ambient = 0.3f};

    const DefaultVertexShader<f32> vertex_shader(model, view, proj);

    renderer.render(mesh.faces, mesh.vertices, mesh.normals, mesh.texture_coordinates, vertex_shader, fragment_shader,
                    depth_buffer, image);

    linear_to_srgb(image);
    dump_ppm(image, std::cout);
}
