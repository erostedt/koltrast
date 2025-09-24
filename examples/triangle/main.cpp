#include <cmath>
#include <cstdlib>
#include <iostream>

#include "camera.hpp"
#include "matrix.hpp"
#include "rasterizer.hpp"
#include "renderer.hpp"
#include "transform.hpp"
#include "types.hpp"

int main()
{
    const std::vector<Vec3f> vertices = {
        {-1.0f, -1.0f, -1.0f}, {2.0f, -1.0f, -3.0f}, {0.0f, 1.0f, -2.0f},
        {1.0f, -1.0f, -1.0f},  {0.0f, 1.0f, -2.0f},  {-2.0f, -1.0f, -3.0f},
    };

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto model = Matrix<f32, 4, 4>::identity();
    const auto view = look_at(Vec3f{0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);
    const auto vp = proj * view;

    std::vector<Vec4f> world_vertices;
    std::vector<Vec3f> world_normals;
    model_to_world(vertices, {}, model, world_vertices, world_normals);

    std::vector<Vec4f> screen_vertices;
    project_to_screen(world_vertices, vp, camera.resolution, screen_vertices);

    ColorImage<f32> image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer<f32>(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);
    rasterize_triangles(screen_vertices, depth_buffer, index_buffer);
    render_triangles(image, {{1.0f, 0, 0}, {0, 1.0f, 0}}, index_buffer);
    dump_ppm(image, std::cout);
    return 0;
}
