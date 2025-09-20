#include <cmath>
#include <cstdlib>
#include <iostream>

#include "camera.hpp"
#include "check.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
#include "transform.hpp"
#include "types.hpp"

std::vector<RGB> random_colors(size_t n)
{
    std::vector<RGB> colors;
    colors.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        colors.push_back({static_cast<u8>(rand() % 256), static_cast<u8>(rand() % 256), static_cast<u8>(rand() % 256)});
    }
    return colors;
}

int main()
{
    const auto mesh = load_obj("/home/eric/git/koltrast/bob/bob_tri.obj");
    const auto colors = random_colors(mesh.faces.size());

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto model = Matrix<f32, 4, 4>::identity();
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);
    rasterize_mesh(mesh, model, view, proj, camera.resolution, depth_buffer, index_buffer);
    draw_triangles(image, colors, index_buffer);
    dump_ppm(image, std::cout);
}
