#pragma once

#include <concepts>
#include <cstdlib>
#include <filesystem>

#include "image.hpp"

namespace fs = std::filesystem;

using Texture = Image<RGB<u8>>;

u8 *load_rgbu8_image(const fs::path &path, size_t &width, size_t &height, size_t &channels, size_t required_components);
f32 *load_rgbf_image(const fs::path &path, size_t &width, size_t &height, size_t &channels, size_t required_components);
void free_image(void *image);

[[nodiscard]] inline Image<RGB<u8>> load_texture(const fs::path &path) noexcept
{
    size_t width, height, channels;
    u8 *data = load_rgbu8_image(path, width, height, channels, 3ul);
    Image<RGB<u8>> texture(static_cast<size_t>(width), static_cast<size_t>(height));

    size_t size = width * height;
    for (size_t i = 0; i < size; ++i)
    {
        auto &e = texture[i];
        e.r = data[3 * i + 0];
        e.g = data[3 * i + 1];
        e.b = data[3 * i + 2];
    }

    free_image(data);
    return texture;
}

template <std::floating_point T> Image<RGB<T>> load_cubemap(const fs::path &path)
{
    size_t width, height, channels;
    f32 *data = load_rgbf_image(path, width, height, channels, 3ul);

    Image<RGB<T>> cubemap((size_t)width, (size_t)height);
    using namespace std;

    for (size_t i = 0; i < (size_t)width * (size_t)height; ++i)
    {
        cubemap[i].r = T(data[3 * i + 0]);
        cubemap[i].g = T(data[3 * i + 1]);
        cubemap[i].b = T(data[3 * i + 2]);
    }
    free_image(data);
    return cubemap;
}
