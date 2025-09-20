#include <cmath>
#include <cstdlib>
#include <iostream>

#include "matrix.hpp"
#include "rasterizer.hpp"
#include "transform.hpp"
#include "types.hpp"

int main()
{
    const std::vector<Vector<f32, 3>> vertices = {
        {-1.0f, -1.0f, -1.0f}, {2.0f, -1.0f, -3.0f}, {0.0f, 1.0f, -2.0f},
        {1.0f, -1.0f, -1.0f},  {0.0f, 1.0f, -2.0f},  {-2.0f, -1.0f, -3.0f},
    };

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto model = Matrix<f32, 4, 4>::identity();
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);
    const auto mvp = proj * view * model;

    std::vector<Vector<f32, 3>> ndc_vertices;
    map_to_ndc(vertices, mvp, camera.resolution, ndc_vertices);

    Triangle t1{ndc_vertices[0], ndc_vertices[1], ndc_vertices[2]};
    Triangle t2{ndc_vertices[3], ndc_vertices[4], ndc_vertices[5]};

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);
    rasterize_triangles({t1, t2}, depth_buffer, index_buffer);
    draw_triangles(image, {{255, 0, 0}, {0, 255, 0}}, index_buffer);
    dump_ppm(image, std::cout);
    return 0;
}
