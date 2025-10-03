#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "camera.hpp"
#include "image.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
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
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);

    const Lights<f32> lights = {
        .ambient = 0.3f,
        .point_lights = {{.position = {0.0f, 1.0f, 2.0f}, .color = {1.0f, 1.0f, 1.0f}, .specular = 0.8f}},
        .directional_lights = {}};

    const DefaultShader<f32> shader = {camera_position, lights, texture, 16.0f};
    const auto vdata = vertex_shader(mesh, model, view, proj, camera.resolution);
    rasterize_triangles(mesh.faces, vdata, depth_buffer, index_buffer);
    render(image, mesh.faces, vdata, mesh.texture_coordinates, shader, index_buffer);

    dump_ppm(image, std::cout);
}
