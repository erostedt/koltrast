#pragma once
#include <concepts>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "check.hpp"
#include "types.hpp"

struct Face
{
    const static size_t INVALID_INDEX = std::numeric_limits<size_t>::max();
    size_t vertex_indices[3]{INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
    size_t texture_indices[3] = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
    size_t normal_indices[3] = {INVALID_INDEX, INVALID_INDEX, INVALID_INDEX};
};

template <std::floating_point T> struct Mesh
{
    std::vector<Vec3<T>> vertices{};
    std::vector<Vec2<T>> texture_coordinates{};
    std::vector<Vec3<T>> normals{};
    std::vector<Face> faces{};
};

template <std::floating_point T> inline Mesh<T> load_obj(const std::filesystem::path &path)
{
    CHECK(std::filesystem::exists(path));
    CHECK(path.extension() == ".obj");

    std::ifstream file(path);
    CHECK(file && file.is_open());
    Mesh<T> mesh{};

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v")
        {
            auto &v = mesh.vertices.emplace_back();
            iss >> v.x() >> v.y() >> v.z();
        }
        else if (prefix == "vn")
        {
            auto &n = mesh.normals.emplace_back();
            iss >> n.x() >> n.y() >> n.z();
        }
        else if (prefix == "vt")
        {
            auto &t = mesh.texture_coordinates.emplace_back();
            iss >> t.x() >> t.y();
        }
        else if (prefix == "f")
        {
            auto &face = mesh.faces.emplace_back();
            std::string f1;
            size_t index = 0;
            while (std::getline(iss, f1, '/'))
            {
                face.vertex_indices[index] = std::stoul(f1) - 1;
                std::getline(iss, f1, '/');
                face.texture_indices[index] = std::stoul(f1) - 1;
                std::getline(iss, f1, ' ');
                face.normal_indices[index] = std::stoul(f1) - 1;
                ++index;
            }
            CHECK(index == 3);
        }
    }
    return mesh;
}
