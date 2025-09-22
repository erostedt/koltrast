#include <cmath>
#include <cstdlib>

#include "texture.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

[[nodiscard]] Image<RGB<f32>> load_texture(const fs::path &path) noexcept
{
    CHECK(fs::exists(path));
    stbi_set_flip_vertically_on_load(1);
    int w, h, c;
    u8 *data = stbi_load(path.c_str(), &w, &h, &c, STBI_rgb);
    CHECK(data != NULL);
    CHECK(c == 3);
    Image<RGB<f32>> texture(static_cast<size_t>(w), static_cast<size_t>(h));

    size_t size = (size_t)w * (size_t)h;
    for (size_t i = 0; i < size; ++i)
    {
        auto &e = texture[i];
        e.r = (f32)data[3 * i + 0] / 255.0f;
        e.g = (f32)data[3 * i + 1] / 255.0f;
        e.b = (f32)data[3 * i + 2] / 255.0f;
    }
    free(data);
    return texture;
}
