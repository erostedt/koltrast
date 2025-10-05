#pragma once

#include <concepts>
#include <execution>
#include <span>

template <typename T> class Buffer3
{
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

  public:
    Buffer3() = default;
    Buffer3(size_t width, size_t height, size_t depth)
        : _width(width), _height(height), _depth(depth), _data(width * height * depth)
    {
    }

    [[nodiscard]] constexpr inline size_t width() const noexcept
    {
        return _width;
    }

    [[nodiscard]] constexpr inline size_t height() const noexcept
    {
        return _height;
    }

    [[nodiscard]] constexpr inline size_t depth() const noexcept
    {
        return _depth;
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return width() * height() * depth();
    }

    [[nodiscard]] constexpr inline T &operator[](size_t i) noexcept
    {
        return _data[i];
    }

    [[nodiscard]] constexpr inline size_t ravel_index(size_t x, size_t y, size_t z) const noexcept
    {
        return (y * width() + x) * depth() + z;
    }

    [[nodiscard]] constexpr inline T &operator[](size_t x, size_t y, size_t z) noexcept
    {
        return operator[](ravel_index(x, y, z));
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t i) const
    {
        return _data[i];
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t x, size_t y, size_t z) const noexcept
    {
        return operator[](ravel_index(x, y, z));
    }

    [[nodiscard]] constexpr inline iterator begin() noexcept
    {
        return std::begin(_data);
    }

    [[nodiscard]] constexpr inline iterator end() noexcept
    {
        return std::end(_data);
    }

    [[nodiscard]] constexpr inline const_iterator begin() const noexcept
    {
        return std::cbegin(_data);
    }

    [[nodiscard]] constexpr inline const_iterator end() const noexcept
    {
        return std::cend(_data);
    }

    [[nodiscard]] constexpr inline std::span<T> slice(size_t x, size_t y) noexcept
    {
        size_t from = ravel_index(x, y, 0);
        size_t to = ravel_index(x, y, depth());
        std::span<T> s(&_data[from], &_data[to]);
        return s;
    }

    [[nodiscard]] constexpr inline std::span<const T> slice(size_t x, size_t y) const noexcept
    {
        size_t from = ravel_index(x, y, 0);
        size_t to = ravel_index(x, y, depth());
        std::span<const T> s(&_data[from], &_data[to]);
        return s;
    }

  private:
    size_t _width;
    size_t _height;
    size_t _depth;
    std::vector<T> _data;
};

template <std::floating_point T> using DepthBuffer = Buffer3<T>;

template <std::floating_point T> constexpr inline void reset_depth_buffer(DepthBuffer<T> &buffer) noexcept
{
    using namespace std;
    fill(execution::par_unseq, begin(buffer), end(buffer), numeric_limits<T>::infinity());
}

template <std::floating_point T>
[[nodiscard]] inline DepthBuffer<T> create_depth_buffer(size_t width, size_t height, size_t depth)
{
    using namespace std;
    DepthBuffer<T> buffer(width, height, depth);
    reset_depth_buffer(buffer);
    return buffer;
}
