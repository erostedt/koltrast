#include <cassert>
#include <cmath>
#include <iostream>

#include "matrix.hpp"
#include "rasterizer.hpp"
#include "transform.hpp"
#include "types.hpp"

int main()
{
    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 1.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    const Vector<f32, 3> a{-1.0f, -1.0f, -1.0f};
    const Vector<f32, 3> b{2.0f, -1.0f, -3.0f};
    const Vector<f32, 3> c{0.0f, 1.0f, -2.0f};

    const Vector<f32, 3> d{1.0f, -1.0f, -1.0f};
    const Vector<f32, 3> e{0.0f, 1.0f, -2.0f};
    const Vector<f32, 3> f{-2.0f, -1.0f, -3.0f};

    const auto as = project_to_screen(a, view, proj, camera.resolution);
    const auto bs = project_to_screen(b, view, proj, camera.resolution);
    const auto cs = project_to_screen(c, view, proj, camera.resolution);

    const auto ds = project_to_screen(d, view, proj, camera.resolution);
    const auto es = project_to_screen(e, view, proj, camera.resolution);
    const auto fs = project_to_screen(f, view, proj, camera.resolution);

    Triangle t1{as, bs, cs};
    Triangle t2{ds, es, fs};

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);
    rasterize_triangles({t1, t2}, depth_buffer, index_buffer);
    draw_triangles(image, {{255, 0, 0}, {0, 255, 0}}, index_buffer);
    dump_ppm(image, std::cout);
}
