#include <X11/Xlib.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "camera.hpp"
#include "counting_iterator.hpp"
#include "depth_buffer.hpp"
#include "image.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "renderer.hpp"
#include "texture.hpp"
#include "transform.hpp"
#include "types.hpp"
#include "window.hpp"

#include <X11/Xutil.h>
#include <cassert>
#include <cstring>
#include <iostream>

using namespace std::chrono;

class PrintFps
{
  public:
    PrintFps(size_t precision = 3) : prec(precision)
    {
        ns = steady_clock::now();
    }

    ~PrintFps()
    {
        const auto diff = steady_clock::now() - ns;
        const auto fps = 1.0 / (1e-9 * (f64)diff.count());
        std::cout << std::setprecision((i32)prec) << fps << std::endl;
    }

  private:
    size_t prec;
    steady_clock::time_point ns;
};

class LimitFps
{
  public:
    LimitFps(f64 target_fps)
    {
        target_elapsed_ns = round<nanoseconds>(duration<f64>(1.0 / target_fps));
        start_ns = steady_clock::now();
    }

    ~LimitFps()
    {
        const auto elapsed_ns = steady_clock::now() - start_ns;
        const auto diff = target_elapsed_ns - elapsed_ns;
        if (diff > 0ns)
        {
            std::this_thread::sleep_for(diff);
        }
    }

  private:
    nanoseconds target_elapsed_ns;
    steady_clock::time_point start_ns;
};

inline void clear_background(Image<RGB<f32>> &image, const Image<RGB<f32>> &cubemap, const Camera<f32> &camera,
                             const Vec3f &camera_position, const Mat4x4f &view) noexcept
{
    // Thid is very slow
    const auto view_port = create_view_port(camera, camera_position, view);
    std::for_each(std::execution::par_unseq, counting_iterator(0), counting_iterator(image.size()), [&](size_t i) {
        size_t x = i % image.width();
        size_t y = i / image.width();

        const Vec3 pixel_center =
            view_port.pixel_top_left + ((f32)x * view_port.pixel_delta_u) + ((f32)y * view_port.pixel_delta_v);
        const Vec3 ray_origin = camera_position;
        const Vec3 ray_direction = *(pixel_center - ray_origin).normalized();
        image[x, y] = sample_cubemap(ray_direction, cubemap);
    });
}

inline void clear_background(Image<RGB<f32>> &image, const RGB<f32> &color) noexcept
{
    using namespace std;
    fill(execution::par_unseq, begin(image), end(image), color);
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << "./mesh obj_file texture_file cubemap_file";
        return 1;
    }

    const auto mesh = load_obj<f32>(argv[1]);
    const auto texture = load_texture(argv[2]);
    const auto cubemap = load_cubemap<f32>(argv[3]);

    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};
    const Vec3f camera_position = {0.0f, 0.5f, 2.0f};
    const auto view = look_at<f32>(camera_position, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    ColorImage<f32> image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer<f32>(camera.resolution.width, camera.resolution.height);

    XWindow window = XWindow::create(camera.resolution.width, camera.resolution.height);

    size_t degrees = 0;
    const DefaultFragmentShader<f32> fragment_shader = {
        .camera_position = camera_position,
        .texture = texture,
        .object_shininess = 16.0f,
        .point_lights = {{.position = {0.0f, 1.0f, 2.0f}, .color = {1.0f, 1.0f, 1.0f}, .specular = 0.8f}},
        .directional_lights = {},
        .ambient = 0.3f};

    auto renderer = create_compatible_renderer(depth_buffer);
    while (!window.should_close)
    {
        PrintFps print_fps;
        // LimitFps limitfps(60.0);
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
        const DefaultVertexShader<f32> vertex_shader(model, view, proj);

        reset_depth_buffer(depth_buffer);
        // clear_background(image, cubemap, camera, camera_position, view);
        clear_background(image, BLACK<f32>);
        DrawFrame frame(window);

        renderer.render(mesh.faces, mesh.vertices, mesh.normals, mesh.texture_coordinates, vertex_shader,
                        fragment_shader, depth_buffer, image);
        linear_to_srgb(image);
        frame.blit(image);
    }
}
