#pragma once
#include <filesystem>

#include "camera.hpp"
#include "types.hpp"

namespace fs = std::filesystem;
using Texture = Image<RGB<f32>>;

[[nodiscard]] Image<RGB<f32>> load_texture(const fs::path &path) noexcept;
