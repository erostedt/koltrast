#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>

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
