#include "camera.hpp"
#include "rasterizer.hpp"
#include "types.hpp"
#include <filesystem>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

RGB<f32> sample_cubemap(Vec3f direction, const Image<RGB<f32>> &cubemap)
{
    f32 phi = std::atan2(direction.z(), direction.x());
    f32 theta = std::acos(direction.y());

    f32 pi = std::numbers::pi_v<f32>;
    f32 u = (phi + pi) / (2.0f * pi);
    f32 v = theta / pi;

    // TODO: (ecrt) Use bilinear interpolation
    size_t px = (size_t)(u * (f32)(cubemap.width() - 1));
    size_t py = (size_t)(v * (f32)(cubemap.height() - 1));

    return cubemap[px, py];
}

Image<RGB<f32>> load_cubemap()
{
    const auto path = fs::path("/home/eric/Downloads/qwantani_sunset_puresky_1k.hdr");
    int width, height, channels;

    stbi_set_flip_vertically_on_load(1);
    assert(fs::exists(path));
    f32 *data = stbi_loadf(path.c_str(), &width, &height, &channels, 0);
    assert(data != NULL);
    assert(channels == 3);

    Image<RGB<f32>> cubemap((size_t)width, (size_t)height);
    using namespace std;

    for (size_t i = 0; i < (size_t)width * (size_t)height; ++i)
    {
        cubemap[i].r = data[3 * i + 0];
        cubemap[i].g = data[3 * i + 1];
        cubemap[i].b = data[3 * i + 2];
    }
    stbi_image_free(data);
    return cubemap;
}

int main()
{
    const auto cubemap = load_cubemap();
    const Camera<f32> camera = {{1280, 720}, 60, 0.2f, 100.0f};

    Image<RGB<f32>> output(1280, 720);
    Vec3f up = {0.0f, 1.0f, 0.0f};
    Vec3f position = {0.0f, 0.0f, 2.0f};
    Vec3f target = {0.0f, 0.0f, 0.0f};

    f32 focus_distance = 1.0f;
    f32 aspect_ratio = (f32)camera.resolution.width / (f32)camera.resolution.height;

    f32 angle = radians(camera.vertical_field_of_view);
    f32 height = std::tan(angle / 2.0f);
    f32 viewport_height = 2.0f * height * focus_distance;
    f32 viewport_width = viewport_height * aspect_ratio;

    Vec3 forward = -*(target - position).normalized();
    Vec3 rightward = *up.cross(forward).normalized();
    Vec3 upward = *forward.cross(rightward).normalized();

    Vec3 viewport_u = viewport_width * rightward;
    Vec3 viewport_v = -viewport_height * upward;

    Vec3 viewport_center = position - focus_distance * forward;
    Vec3 viewport_upper_left = viewport_center - viewport_u / 2.0f - viewport_v / 2.0f;
    Vec3 PixelDeltaU = viewport_u / (f32)camera.resolution.width;
    Vec3 PixelDeltaV = viewport_v / (f32)camera.resolution.height;
    Vec3 PixelTopLeft = viewport_upper_left + 0.5f * (PixelDeltaU + PixelDeltaV);

    for (size_t y = 0; y < camera.resolution.height; ++y)
    {
        for (size_t x = 0; x < camera.resolution.width; ++x)
        {
            const f32 fx = (f32)x + 0.5f;
            const f32 fy = (f32)y + 0.5f;
            const Vec3 pixel_center = PixelTopLeft + (fx * PixelDeltaU) + (fy * PixelDeltaV);
            const Vec3 ray_origin = position;
            const Vec3 ray_direction = *(pixel_center - ray_origin).normalized();
            output[x, y] = sample_cubemap(ray_direction, cubemap);
        }
    }
    dump_ppm(output, std::cout);
}
