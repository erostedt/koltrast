#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>

template <typename T, typename U, size_t N, typename BinaryOperation>
constexpr auto Elementwise(const std::array<T, N> &left, const std::array<U, N> &right, BinaryOperation op)
{
    using namespace std;
    using value_type = decltype(op(std::declval<T>(), std::declval<U>()));
    std::array<value_type, N> out;
    transform(begin(left), end(left), begin(right), begin(out), op);
    return out;
}

template <typename T, typename U, size_t N, typename BinaryOperation>
constexpr auto Elementwise(const std::array<T, N> &left, U right, BinaryOperation op)
{
    using namespace std;
    using value_type = decltype(op(std::declval<T>(), std::declval<U>()));
    std::array<value_type, N> out;
    transform(begin(left), end(left), begin(out), [&right, &op](const T &element) { return op(element, right); });
    return out;
}

template <typename T, typename U, size_t N, typename BinaryOperation>
constexpr auto Elementwise(T left, const std::array<U, N> &right, BinaryOperation op)
{
    using namespace std;
    using value_type = decltype(op(std::declval<T>(), std::declval<U>()));
    std::array<value_type, N> out;
    transform(begin(right), end(right), begin(out), [&left, &op](const U &element) { return op(left, element); });
    return out;
}

template <typename T, size_t N> class _Vector
{
    using value_type = std::array<T, N>::value_type;
    using pointer = std::array<T, N>::pointer;
    using const_pointer = std::array<T, N>::const_pointer;
    using reference = std::array<T, N>::reference;
    using const_reference = std::array<T, N>::const_reference;
    using iterator = std::array<T, N>::iterator;
    using const_iterator = std::array<T, N>::const_iterator;
    using size_type = std::array<T, N>::size_type;
    using difference_type = std::array<T, N>::difference_type;
    using reverse_iterator = std::array<T, N>::reverse_iterator;
    using const_reverse_iterator = std::array<T, N>::const_reverse_iterator;

  public:
    template <typename... U> constexpr _Vector(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == N, "Wrong number of elements");
        static_assert((std::is_same_v<U, T> && ...), "All arguments must be exactly T");
    }

    constexpr _Vector(const std::array<T, N> &arr) : elements(arr)
    {
    }

    [[nodiscard]] constexpr inline iterator begin() noexcept
    {
        return std::begin(elements);
    }

    [[nodiscard]] constexpr inline const_iterator begin() const noexcept
    {
        return std::begin(elements);
    }

    [[nodiscard]] constexpr inline iterator end() noexcept
    {
        return std::end(elements);
    }

    [[nodiscard]] constexpr inline const_iterator end() const noexcept
    {
        return std::end(elements);
    }

    [[nodiscard]] constexpr inline bool operator==(const _Vector &v) const noexcept
    {
        return std::equal(begin(), end(), std::begin(v));
    }

    [[nodiscard]] constexpr inline bool operator!=(const _Vector &v) const noexcept
    {
        return !(*this == v);
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return N;
    }

    [[nodiscard]] constexpr inline T dot(const _Vector &v) const noexcept
    {
        return std::inner_product(this->begin(), this->end(), std::begin(v), T{});
    }

    [[nodiscard]] constexpr inline T squared_length() const noexcept
    {
        return this->dot(*this);
    }

    [[nodiscard]] constexpr inline T length() const noexcept
    {
        return std::sqrt(squared_length());
    }

    std::array<T, N> elements{};
};

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator+(const _Vector<T, N> &left, const _Vector<T, N> &right) noexcept
{
    return Elementwise(left.elements, right.elements, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator+(const _Vector<T, N> &left, const T &right) noexcept
{
    return Elementwise(left.elements, right, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator+(const T &left, const _Vector<T, N> &right) noexcept
{
    return Elementwise(left, right.elements, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator-(const _Vector<T, N> &left, const _Vector<T, N> &right) noexcept
{
    return Elementwise(left.elements, right.elements, std::minus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator-(const _Vector<T, N> &left, const T &right) noexcept
{
    return Elementwise(left.elements, right, std::minus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator-(const T &left, const _Vector<T, N> &right) noexcept
{
    return Elementwise(left, right.elements, std::minus<T>{});
}

template <typename T, size_t N> [[nodiscard]] constexpr inline _Vector<T, N> operator-(const _Vector<T, N> &v) noexcept
{
    using namespace std;
    auto out = v;
    transform(begin(v), end(v), begin(out), std::negate<T>{});
    return out;
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator*(const _Vector<T, N> &left, T right) noexcept
{
    return Elementwise(left.elements, right, std::multiplies<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator*(T left, const _Vector<T, N> &right) noexcept
{
    return Elementwise(left, right.elements, std::multiplies<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline _Vector<T, N> operator/(const _Vector<T, N> &numerator, T denominator) noexcept
{
    return Elementwise(numerator.elements, denominator, std::divides<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline T squared_distance(const _Vector<T, N> &from, const _Vector<T, N> &to) noexcept
{
    return (from - to).squared_length();
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline T distance(const _Vector<T, N> &from, const _Vector<T, N> &to) noexcept
{
    return std::sqrt(squared_distance(from, to));
}

template <typename T, std::size_t N> class Vector : public _Vector<T, N>
{
  public:
    using _Vector<T, N>::_Vector;
};

template <typename T> class Vector<T, 2> : public _Vector<T, 2>
{
  public:
    using _Vector<T, 2>::_Vector;

    [[nodiscard]] constexpr inline T x() const noexcept
    {
        return this->elements[0];
    }

    [[nodiscard]] constexpr inline T y() const noexcept
    {
        return this->elements[1];
    }
};

template <typename T> class Vector<T, 3> : public _Vector<T, 3>
{
  public:
    using _Vector<T, 3>::_Vector;

    [[nodiscard]] constexpr inline T x() const noexcept
    {
        return this->elements[0];
    }

    [[nodiscard]] constexpr inline T y() const noexcept
    {
        return this->elements[1];
    }

    [[nodiscard]] constexpr inline T z() const noexcept
    {
        return this->elements[2];
    }

    [[nodiscard]] constexpr inline Vector cross(const Vector &v) const noexcept
    {
        const auto &u = *this;
        const auto x = u.y() * v.z() - u.z() * v.y();
        const auto y = u.z() * v.x() - u.x() * v.z();
        const auto z = u.x() * v.y() - u.y() * v.x();
        return {x, y, z};
    }
};

template <typename T> class Vector<T, 4> : public _Vector<T, 4>
{
  public:
    using _Vector<T, 4>::_Vector;

    [[nodiscard]] constexpr inline T x() const noexcept
    {
        return this->elements[0];
    }

    [[nodiscard]] constexpr inline T y() const noexcept
    {
        return this->elements[1];
    }

    [[nodiscard]] constexpr inline T z() const noexcept
    {
        return this->elements[2];
    }

    [[nodiscard]] constexpr inline T w() const noexcept
    {
        return this->elements[3];
    }
};
