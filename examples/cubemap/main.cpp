#include <filesystem>

#include "camera.hpp"
#include "shader.hpp"
#include "texture.hpp"
#include "transform.hpp"
#include "types.hpp"

namespace fs = std::filesystem;

int main()
{
    // TODO: (ecrt) Remove hardcoded path
    const auto path = fs::path("/home/eric/Downloads/qwantani_sunset_puresky_1k.hdr");
    const auto cubemap = load_cubemap<f32>(path);
    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};

    Vec3f up = {0.0f, 1.0f, 0.0f};
    Vec3f position = {0.0f, 0.0f, 2.0f};
    Vec3f target = {0.0f, 0.0f, 0.0f};
    const auto view = look_at<f32>(position, target, up);
    const auto v = create_view_port(camera, position, view);

    Image<RGB<f32>> output(1280, 720);

    for (size_t y = 0; y < camera.resolution.height; ++y)
    {
        for (size_t x = 0; x < camera.resolution.width; ++x)
        {
            const f32 fx = (f32)x + 0.5f;
            const f32 fy = (f32)y + 0.5f;
            const Vec3 pixel_center = v.pixel_top_left + (fx * v.pixel_delta_u) + (fy * v.pixel_delta_v);
            const Vec3 ray_origin = position;
            const Vec3 ray_direction = *(pixel_center - ray_origin).normalized();
            output[x, y] = sample_cubemap(ray_direction, cubemap);
        }
    }
    dump_ppm(output, std::cout);
}
