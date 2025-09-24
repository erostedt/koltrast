#include <X11/Xlib.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "camera.hpp"
#include "counting_iterator.hpp"
#include "light.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
#include "renderer.hpp"
#include "texture.hpp"
#include "transform.hpp"
#include "types.hpp"
#include "window.hpp"

#include <X11/Xutil.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

class PrintFps
{

  public:
    PrintFps()
    {
        ns = std::chrono::steady_clock::now();
    }

    ~PrintFps()
    {
        const auto diff = std::chrono::steady_clock::now() - ns;
        const auto fps = 1.0 / (1e-9 * (f64)diff.count());
        std::cout << std::setprecision(3) << fps << std::endl;
    }

  private:
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> ns;
};

void apply_lighting(ColorImage<f32> &image, const Vec3f &camera_position, const std::vector<Face> &faces,
                    const std::vector<Vec4f> &screen_vertices, const std::vector<Vec4f> &world_vertices,
                    const std::vector<Vec3f> &world_normals, const IndexBuffer &index_buffer)
{
    PointLight point_light{
        .position = {1.0f, 0.0f, 0.0f},
        .color = {0.0f, 1.0f, 0.0f},
        .specular = 0.5f,
    };

    f32 ambient = 0.1f;
    f32 shininess = 8.0f;
    using namespace std;
    for_each(execution::par_unseq, counting_iterator(0), counting_iterator(index_buffer.size()), [&](size_t i) {
        size_t x = i % index_buffer.width();
        size_t y = i / index_buffer.width();
        size_t index = index_buffer[x, y];
        if (index != std::numeric_limits<size_t>::max())
        {
            const auto &face = faces[index];
            const auto sv1 = screen_vertices[face.vertex_indices[0]];
            const auto sv2 = screen_vertices[face.vertex_indices[1]];
            const auto sv3 = screen_vertices[face.vertex_indices[2]];

            const auto wv1 = world_vertices[face.vertex_indices[0]];
            const auto wv2 = world_vertices[face.vertex_indices[1]];
            const auto wv3 = world_vertices[face.vertex_indices[2]];

            const auto wn1 = world_normals[face.normal_indices[0]];
            const auto wn2 = world_normals[face.normal_indices[1]];
            const auto wn3 = world_normals[face.normal_indices[2]];

            const auto bary = barycentric({(f32)x + 0.5f, (f32)y + 0.5f}, sv1, sv2, sv3);
            const auto world_position = (bary.x() * wv1 + bary.y() * wv2 + bary.z() * wv3).xyz();
            const auto world_normal = (bary.x() * wn1 + bary.y() * wn2 + bary.z() * wn3);

            const RGB<f32> light = sample_light(world_position, world_normal, camera_position, shininess, point_light);

            const RGB<f32> object_color = image[x, y];
            image[x, y] = {(ambient + light.r) * object_color.r, (ambient + light.g) * object_color.g,
                           (ambient + light.b) * object_color.b};
        }
    });
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << "./mesh obj_file texture_file > output.ppm";
        return 1;
    }

    const auto mesh = load_obj<f32>(argv[1]);
    const auto texture = load_texture<f32>(argv[2]);

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const Vec3f camera_position = {0.0f, 0.0f, 2.0f};
    const auto view = look_at<f32>(camera_position, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    ColorImage<f32> image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer<f32>(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);

    XWindow window = XWindow::create(camera.resolution.width, camera.resolution.height);

    std::vector<Vec4f> world_vertices;
    std::vector<Vec3f> world_normals;
    std::vector<Vec4f> screen_vertices;
    size_t degrees = 0;
    while (!window.should_close)
    {
        PrintFps fps;
        auto events = window.poll_events();
        for (auto &event : events)
        {
            if (event.type == KeyPress && XLookupKeysym(&event.xkey, 0) == 'q')
            {
                window.should_close = true;
            }
        }

        degrees = (degrees + 1) % 360;
        const auto model = model_matrix<f32>({0.0f, 0.0f, 0.0f}, {0.0f, (f32)degrees, 0.0f}, {1.0f, 1.0f, 1.0f});
        const auto vp = proj * view;

        reset_depth_buffer(depth_buffer);
        reset_index_buffer(index_buffer);
        std::fill(image.begin(), image.end(), BLACK<f32>);

        DrawFrame frame(window);

        model_to_world(mesh.vertices, mesh.normals, model, world_vertices, world_normals);

        project_to_screen(world_vertices, vp, camera.resolution, screen_vertices);
        rasterize_triangles(mesh.faces, screen_vertices, depth_buffer, index_buffer);
        render_triangles(image, mesh.faces, screen_vertices, mesh.texture_coordinates, texture, index_buffer);
        apply_lighting(image, camera_position, mesh.faces, screen_vertices, world_vertices, world_normals,
                       index_buffer);

        frame.blit(image);
    }
}
