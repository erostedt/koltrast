#pragma once
#include "camera.hpp"
#include "matrix.hpp"
#include "obj.hpp"
#include "types.hpp"

using ColorImage = Image<RGB>;
using DepthBuffer = Image<f32>;
using IndexBuffer = Image<size_t>;

struct Triangle
{
    Vector<f32, 3> p1;
    Vector<f32, 3> p2;
    Vector<f32, 3> p3;
};

struct BoundingBox
{
    Vector<size_t, 2> top_left;
    Vector<size_t, 2> bottom_right;
};

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height);
[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height);

void rasterize_triangles(const std::vector<Triangle> &screen_triangles, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept;
void rasterize_triangles(const std::vector<Face> &faces, const std::vector<Vector<f32, 4>> &screen_vertices,
                         DepthBuffer &depth_buffer, IndexBuffer &index_buffer) noexcept;

void draw_triangles(ColorImage &image, const std::vector<RGB> &colors, const IndexBuffer &index_buffer) noexcept;
void draw_triangles(ColorImage &image, const std::vector<Face> &faces,
                    const std::vector<Vector<f32, 4>> &screen_vertices,
                    const std::vector<Vector<f32, 2>> &texture_coordinates, const Image<RGB> &texture,
                    const IndexBuffer &index_buffer) noexcept;
void dump_ppm(const ColorImage &image, std::ostream &stream);
