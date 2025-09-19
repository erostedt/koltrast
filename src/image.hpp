#pragma once

#include <cmath>
#include <vector>

#include "types.hpp"

struct RGB
{
    u8 r;
    u8 g;
    u8 b;
};

template <typename T> class Image
{
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

  public:
    Image(u32 width, u32 height) : m_width(width), m_height(height), pixels(width * height)
    {
    }

    inline u32 width() const
    {
        return m_width;
    }

    inline u32 height() const
    {
        return m_height;
    }

    inline T &operator[](u32 i)
    {
        return pixels[i];
    }

    inline T &operator[](u32 x, u32 y)
    {
        return operator[](y * m_width + x);
    }

    inline const T &operator[](u32 i) const
    {
        return pixels[i];
    }

    inline const T &operator[](u32 x, u32 y) const
    {
        return operator[](y * m_width + x);
    }

    inline iterator begin()
    {
        return std::begin(pixels);
    }

    inline iterator end()
    {
        return std::end(pixels);
    }

    inline const_iterator begin() const
    {
        return std::cbegin(pixels);
    }

    inline const_iterator end() const
    {
        return std::cend(pixels);
    }

  private:
    u32 m_width;
    u32 m_height;
    std::vector<T> pixels;
};
