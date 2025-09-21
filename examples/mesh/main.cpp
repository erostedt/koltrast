#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "camera.hpp"
#include "check.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
#include "transform.hpp"
#include "types.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

Image<RGB<f32>> load_texture(const fs::path &path)
{
    CHECK(fs::exists(path));
    stbi_set_flip_vertically_on_load(1);
    int w, h, c;
    u8 *data = stbi_load(path.c_str(), &w, &h, &c, STBI_rgb);
    CHECK(data != NULL);
    CHECK(c == 3);
    Image<RGB<f32>> texture(static_cast<size_t>(w), static_cast<size_t>(h));

    size_t size = (size_t)w * (size_t)h;
    for (size_t i = 0; i < size; ++i)
    {
        auto &e = texture[i];
        e.r = (f32)data[3 * i + 0] / 255.0f;
        e.g = (f32)data[3 * i + 1] / 255.0f;
        e.b = (f32)data[3 * i + 2] / 255.0f;
    }
    free(data);
    return texture;
}

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
