#include <cmath>
#include <cstdlib>
#include <filesystem>

#include "check.hpp"
#include "types.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace fs = std::filesystem;

u8 *load_image(const fs::path &path, size_t &width, size_t &height, size_t &channels, size_t required_components)
{
    CHECK(fs::exists(path));
    stbi_set_flip_vertically_on_load(1);
    i32 w, h, c;
    u8 *data = stbi_load(path.c_str(), &w, &h, &c, (i32)required_components);
    CHECK(data != NULL);
    CHECK(c == 3);
    width = (size_t)w;
    height = (size_t)h;
    channels = (size_t)c;
    return data;
}
