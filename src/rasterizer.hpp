#pragma once
#include "camera.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "texture.hpp"
#include "types.hpp"

using ColorImage = Image<RGB<u8>>;
using DepthBuffer = Image<f32>;
using IndexBuffer = Image<size_t>;

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height);
void reset_depth_buffer(DepthBuffer &buffer) noexcept;

[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height);
void reset_index_buffer(IndexBuffer &buffer) noexcept;

void rasterize_triangles(const std::vector<Vector<f32, 4>> &screen_vertices, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept;
void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vector<f32, 4>> &screen_vertices,
                         DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept;

void draw_triangles(ColorImage &image, const std::vector<RGB<u8>> &colors, const IndexBuffer &index_buffer) noexcept;
void draw_triangles(ColorImage &image, const std::vector<Face> &faces,
                    const std::vector<Vector<f32, 4>> &screen_vertices,
                    const std::vector<Vector<f32, 2>> &texture_coordinates, const Texture &texture,
                    const IndexBuffer &index_buffer) noexcept;
void dump_ppm(const ColorImage &image, std::ostream &stream);
