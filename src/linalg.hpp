#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <numeric>

template <typename T, typename U>
concept Addable = requires(T a, U b) {
    { a + b };
};

template <typename T, typename U>
concept Subtractable = requires(T a, U b) {
    { a - b };
};

template <typename T>
concept Negatable = requires(T a) {
    { -a } -> std::same_as<T>;
};

template <typename T, typename U>
concept Multiplicable = requires(T a, U b) {
    { a * b };
};

template <typename T, typename U>
concept Divisable = requires(T a, U b) {
    { a / b };
};

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

template <typename T, size_t N>
    requires std::is_signed_v<T>
class Vector
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
    template <typename... U> constexpr Vector(U... elems) : elements{elems...}
    {
        static_assert(sizeof...(U) == N, "Wrong number of elements");
        static_assert((std::is_same_v<U, T> && ...), "All arguments must be exactly T");
    }

    constexpr Vector(const std::array<T, N> &arr) : elements(arr)
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

    [[nodiscard]] constexpr inline bool operator==(const Vector &v) const noexcept
        requires std::equality_comparable<T>
    {
        return elements == v.elements;
    }

    [[nodiscard]] constexpr inline bool operator!=(const Vector &v) const noexcept
        requires std::equality_comparable<T>
    {
        return !(*this == v);
    }

    [[nodiscard]] constexpr inline size_t size() const noexcept
    {
        return N;
    }

    [[nodiscard]] constexpr inline T dot(const Vector &v) const noexcept
        requires Multiplicable<T, T> && Addable<T, T>
    {
        return std::inner_product(this->begin(), this->end(), std::begin(v), T{});
    }

    [[nodiscard]] constexpr inline T squared_length() const noexcept
        requires std::floating_point<T>
    {
        return this->dot(*this);
    }

    [[nodiscard]] constexpr inline T length() const noexcept
        requires std::floating_point<T>
    {
        return std::sqrt(squared_length());
    }

    [[nodiscard]] constexpr inline T x() const noexcept
        requires(N >= 1 && N <= 4)
    {
        return elements[0];
    }

    [[nodiscard]] constexpr inline T y() const noexcept
        requires(N >= 2 && N <= 4)
    {
        return elements[1];
    }

    [[nodiscard]] constexpr inline T z() const noexcept
        requires(N >= 3 && N <= 4)
    {
        return elements[2];
    }

    [[nodiscard]] constexpr inline T w() const noexcept
        requires(N == 4)
    {
        return elements[3];
    }

    [[nodiscard]] constexpr inline Vector cross(const Vector &v) const noexcept
        requires(N == 3) && Multiplicable<T, T> && Subtractable<T, T>
    {
        const auto &u = *this;
        const auto x = u.y() * v.z() - u.z() * v.y();
        const auto y = u.z() * v.x() - u.x() * v.z();
        const auto z = u.x() * v.y() - u.y() * v.x();
        return {x, y, z};
    }

    std::array<T, N> elements{};
};

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator+(const Vector<T, N> &left, const Vector<T, N> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right.elements, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator+(const Vector<T, N> &left, const T &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left.elements, right, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator+(const T &left, const Vector<T, N> &right) noexcept
    requires Addable<T, T>
{
    return Elementwise(left, right.elements, std::plus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator-(const Vector<T, N> &left, const Vector<T, N> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right.elements, std::minus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator-(const Vector<T, N> &left, const T &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left.elements, right, std::minus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator-(const T &left, const Vector<T, N> &right) noexcept
    requires Subtractable<T, T>
{
    return Elementwise(left, right.elements, std::minus<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator-(const Vector<T, N> &v) noexcept
    requires Negatable<T>
{
    using namespace std;
    auto out = v;
    transform(begin(v), end(v), begin(out), std::negate<T>{});
    return out;
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator*(const Vector<T, N> &left, T right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left.elements, right, std::multiplies<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator*(T left, const Vector<T, N> &right) noexcept
    requires Multiplicable<T, T>
{
    return Elementwise(left, right.elements, std::multiplies<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline Vector<T, N> operator/(const Vector<T, N> &numerator, T denominator) noexcept
    requires Divisable<T, T>
{
    return Elementwise(numerator.elements, denominator, std::divides<T>{});
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline T squared_distance(const Vector<T, N> &from, const Vector<T, N> &to) noexcept
    requires std::floating_point<T>
{
    return (from - to).squared_length();
}

template <typename T, size_t N>
[[nodiscard]] constexpr inline T distance(const Vector<T, N> &from, const Vector<T, N> &to) noexcept
    requires std::floating_point<T>
{
    return std::sqrt(squared_distance(from, to));
}
