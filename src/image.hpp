#pragma once

#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <vector>

#include "color.hpp"

struct Resolution
{
    size_t width = 640;
    size_t height = 480;
};

inline bool operator==(const Resolution &lhs, const Resolution &rhs)
{
    return (lhs.width == rhs.width) && (lhs.height == rhs.height);
}

inline bool operator!=(const Resolution &lhs, const Resolution &rhs)
{
    return !(lhs == rhs);
}

template <typename T> class Image
{
    using iterator = std::vector<T>::iterator;
    using const_iterator = std::vector<T>::const_iterator;

  public:
    Image() = default;
    Image(size_t width, size_t height) : res(width, height), pixels(width * height)
    {
    }

    [[nodiscard]] constexpr inline size_t width() const noexcept
    {
        return res.width;
    }

    [[nodiscard]] constexpr inline size_t height() const noexcept
    {
        return res.height;
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return width() * height();
    }

    [[nodiscard]] constexpr inline const Resolution &resolution() const noexcept
    {
        return res;
    }

    [[nodiscard]] constexpr inline T &operator[](size_t i) noexcept
    {
        return pixels[i];
    }

    [[nodiscard]] constexpr inline T &operator[](size_t x, size_t y) noexcept
    {
        return operator[](y * width() + x);
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t i) const
    {
        return pixels[i];
    }

    [[nodiscard]] constexpr inline const T &operator[](size_t x, size_t y) const noexcept
    {
        return operator[](y * width() + x);
    }

    [[nodiscard]] constexpr inline iterator begin() noexcept
    {
        return std::begin(pixels);
    }

    [[nodiscard]] constexpr inline iterator end() noexcept
    {
        return std::end(pixels);
    }

    [[nodiscard]] constexpr inline const_iterator begin() const noexcept
    {
        return std::cbegin(pixels);
    }

    [[nodiscard]] constexpr inline const_iterator end() const noexcept
    {
        return std::cend(pixels);
    }

  private:
    Resolution res;
    std::vector<T> pixels;
};

template <std::floating_point T> using ColorImage = Image<RGB<T>>;
template <std::floating_point T> inline void linear_to_srgb(ColorImage<T> &linear) noexcept
{
    using namespace std;
    std::transform(execution::par, begin(linear), end(linear), begin(linear),
                   [](const RGB<T> &color) -> RGB<T> { return linear_to_srgbf(color); });
}

template <std::floating_point T> inline void dump_ppm(const ColorImage<T> &image, std::ostream &stream)
{
    using std::ranges::for_each;
    stream << "P3\n";
    stream << image.width() << ' ' << image.height() << "\n255\n";

    const auto write_pixel = [&stream](const auto &color) {
        const auto r = std::clamp(color.r * T{255}, T{0}, T{255});
        const auto g = std::clamp(color.g * T{255}, T{0}, T{255});
        const auto b = std::clamp(color.b * T{255}, T{0}, T{255});
        stream << (int)r << ' ' << (int)g << ' ' << (int)b << '\n';
    };

    for_each(image, write_pixel);
}
