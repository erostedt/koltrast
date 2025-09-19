#include <cmath>
#include <cstdlib>
#include <iostream>

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

std::vector<Triangle> make_triangles(const Mesh &m)
{
    std::vector<Triangle> tris;
    tris.reserve(m.faces.size());
    for (const auto &face : m.faces)
    {
        tris.push_back({
            m.vertices[face.vertex_indices[0]],
            m.vertices[face.vertex_indices[1]],
            m.vertices[face.vertex_indices[2]],
        });
    }
    return tris;
}

std::vector<Triangle> project_triangles(const std::vector<Triangle> &tris, const Camera<f32> &camera,
                                        const Mat4x4<f32> &view)
{
    const auto proj = projection_matrix(camera);
    std::vector<Triangle> out{};
    out.reserve(tris.size());

    for (const auto &tri : tris)
    {
        const auto a = project_to_screen(tri.p1, view, proj, camera.resolution);
        const auto b = project_to_screen(tri.p2, view, proj, camera.resolution);
        const auto c = project_to_screen(tri.p3, view, proj, camera.resolution);
        out.push_back({a, b, c});
    }
    return out;
}

int main()
{
    const auto mesh = load_obj("/home/eric/git/koltrast/bob/bob_tri.obj");
    const auto tris = make_triangles(mesh);
    const auto colors = random_colors(mesh.faces.size());

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const auto view = look_at(Vector<f32, 3>{0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});

    const auto triangles = project_triangles(tris, camera, view);

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);
    rasterize_triangles(triangles, depth_buffer, index_buffer);
    draw_triangles(image, colors, index_buffer);
    dump_ppm(image, std::cout);
}
