#pragma once
#include "camera.hpp"
#include "obj.hpp"
#include "texture.hpp"
#include "types.hpp"

using ColorImage = Image<RGB<f32>>;
using DepthBuffer = Image<f32>;
using IndexBuffer = Image<size_t>;

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height);
void reset_depth_buffer(DepthBuffer &buffer) noexcept;

[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height);
void reset_index_buffer(IndexBuffer &buffer) noexcept;

void rasterize_triangles(const std::vector<Vec4f> &screen_vertices, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept;
void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                         DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept;

void draw_triangles(ColorImage &image, const std::vector<RGB<f32>> &colors, const IndexBuffer &index_buffer) noexcept;
void draw_triangles(ColorImage &image, const std::vector<Face> &faces, const std::vector<Vec4f> &screen_vertices,
                    const std::vector<Vec2f> &texture_coordinates, const Texture &texture,
                    const IndexBuffer &index_buffer) noexcept;
void dump_ppm(const ColorImage &image, std::ostream &stream);
