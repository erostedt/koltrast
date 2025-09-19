#pragma once
#include "camera.hpp"
#include "matrix.hpp"
#include "types.hpp"

using Vec3f = Vector<f32, 3>;
using Vec2s = Vector<size_t, 2>;
using ColorImage = Image<RGB>;
using DepthBuffer = Image<f32>;
using IndexBuffer = Image<size_t>;

struct Triangle
{
    Vec3f p1;
    Vec3f p2;
    Vec3f p3;
};

struct BoundingBox
{
    Vec2s top_left;
    Vec2s bottom_right;
};

[[nodiscard]] DepthBuffer create_depth_buffer(size_t width, size_t height);
[[nodiscard]] IndexBuffer create_index_buffer(size_t width, size_t height);
void rasterize_triangles(const std::vector<Triangle> &triangles, DepthBuffer &depth_buffer,
                         IndexBuffer &index_buffer) noexcept;
void draw_triangles(ColorImage &image, const std::vector<RGB> &colors, const IndexBuffer &index_buffer) noexcept;
void dump_ppm(const ColorImage &image, std::ostream &stream);
