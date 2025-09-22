#include <X11/Xlib.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "camera.hpp"
#include "check.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "rasterizer.hpp"
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
    const auto view = look_at<f32>({0.0f, 0.0f, 2.0f}, {0.0f, 0.0f, -1.0f}, {0.0f, 1.0f, 0.0f});
    const auto proj = projection_matrix(camera);

    ColorImage image(camera.resolution.width, camera.resolution.height);
    auto depth_buffer = create_depth_buffer(camera.resolution.width, camera.resolution.height);
    auto index_buffer = create_index_buffer(camera.resolution.width, camera.resolution.height);

    XWindow window = XWindow::create(camera.resolution.width, camera.resolution.height);
    std::vector<Vector<f32, 4>> screen_vertices;

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
        const auto mvp = proj * view * model;

        reset_depth_buffer(depth_buffer);
        reset_index_buffer(index_buffer);
        std::fill(image.begin(), image.end(), BLACK<u8>);

        DrawFrame frame(window);

        project_to_screen(mesh.vertices, mvp, camera.resolution, screen_vertices);
        rasterize_triangles(mesh.faces, screen_vertices, depth_buffer, index_buffer);
        draw_triangles(image, mesh.faces, screen_vertices, mesh.texture_coordinates, texture, index_buffer);

        frame.blit(image);
    }
}
