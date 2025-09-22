#include <cstddef>
#include <iterator>

struct counting_iterator
{
    using value_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::random_access_iterator_tag;
    using reference = value_type;
    using pointer = void;

    value_type value = 0;

    // Core ops
    value_type operator*() const noexcept
    {
        return value;
    }

    counting_iterator &operator++() noexcept
    {
        ++value;
        return *this;
    }
    counting_iterator operator++(int) noexcept
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    counting_iterator &operator--() noexcept
    {
        --value;
        return *this;
    }
    counting_iterator operator--(int) noexcept
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    counting_iterator &operator+=(difference_type n) noexcept
    {
        value = static_cast<value_type>(static_cast<difference_type>(value) + n);
        return *this;
    }
    counting_iterator &operator-=(difference_type n) noexcept
    {
        value = static_cast<value_type>(static_cast<difference_type>(value) - n);
        return *this;
    }

    friend counting_iterator operator+(counting_iterator it, difference_type n) noexcept
    {
        return it += n;
    }
    friend counting_iterator operator+(difference_type n, counting_iterator it) noexcept
    {
        return it += n;
    }
    friend counting_iterator operator-(counting_iterator it, difference_type n) noexcept
    {
        return it -= n;
    }
    friend difference_type operator-(counting_iterator a, counting_iterator b) noexcept
    {
        return static_cast<difference_type>(a.value) - static_cast<difference_type>(b.value);
    }

    value_type operator[](difference_type n) const noexcept
    {
        return static_cast<value_type>(static_cast<difference_type>(value) + n);
    }

    // C++20 comparisons
    auto operator<=>(const counting_iterator &) const = default;
};
