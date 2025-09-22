#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "camera.hpp"
#include "check.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
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

    const auto mesh = load_obj(argv[1]);
    const auto texture = load_texture(argv[2]);

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto model = model_matrix<f32>({0.0f, 0.0f, 0.0f}, {0.0f, 45.0f, 0.0f}, {1.0f, 1.0f, 1.0f});
    const auto view = look_at<f32>({0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);
    const auto mvp = proj * view * model;

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);

    std::vector<Vector<f32, 4>> screen_vertices;
    project_to_screen(mesh.vertices, mvp, camera.resolution, screen_vertices);
    rasterize_triangles(mesh.faces, screen_vertices, depth_buffer, index_buffer);
    draw_triangles(image, mesh.faces, screen_vertices, mesh.texture_coordinates, texture, index_buffer);

    dump_ppm(image, std::cout);
}
