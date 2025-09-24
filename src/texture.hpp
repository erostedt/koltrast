#pragma once
#include <concepts>
#include <filesystem>

#include "camera.hpp"

namespace fs = std::filesystem;

template <std::floating_point T> using Texture = Image<RGB<T>>;

u8 *load_image(const fs::path &path, size_t &width, size_t &height, size_t &channels, size_t required_components);

template <std::floating_point T> [[nodiscard]] Image<RGB<T>> load_texture(const fs::path &path) noexcept
{
    size_t width, height, channels;
    u8 *data = load_image(path, width, height, channels, 3ul);
    Image<RGB<T>> texture(static_cast<size_t>(width), static_cast<size_t>(height));

    size_t size = width * height;
    for (size_t i = 0; i < size; ++i)
    {
        auto &e = texture[i];
        e.r = (T)data[3 * i + 0] / T{255};
        e.g = (T)data[3 * i + 1] / T{255};
        e.b = (T)data[3 * i + 2] / T{255};
    }
    free(data);
    return texture;
}
