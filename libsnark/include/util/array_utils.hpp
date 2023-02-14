#pragma once

#include <array>
#include <cstddef>

template<typename T, size_t n, typename... Args>
auto make_uniform_array(Args &&...args)
{
    return []<size_t... i>(Args && ...args, std::index_sequence<i...>)
    {
        return std::array{((void)i, T{std::forward<Args>(args)...})...};
    }
    (std::forward<Args>(args)..., std::make_index_sequence<n>{});
}

template<typename T, typename... Args>
auto make_uniform_array(Args &&...args)
{
    return make_uniform_array<typename T::value_type, std::tuple_size_v<T>>(
        std::forward<Args>(args)...);
}

template<typename T, size_t n>
struct ArrayPP : public std::array<T, n>
{
    /* ArrayPP
    * This class extends the standard array by allowing (compile-time) variadic size initialization
    * from some arguments. 
    * This is particularly useful if move-construction of T is deleted or expensive, as no move 
    * or copy is required. 
    * If moving is not a concern, then the make_uniform_array function will work just fine. 
    */
    using super = std::array<T, n>;
    using super::const_iterator;
    using super::const_pointer;
    using super::const_reference;
    using super::const_reverse_iterator;
    using super::difference_type;
    using super::iterator;
    using super::pointer;
    using super::reference;
    using super::reverse_iterator;
    using super::size_type;
    using super::value_type;


    template<typename... Args>
    static auto init(Args &&...args)
    {
        if constexpr (sizeof...(Args) < n && std::is_constructible_v<T, Args...>)
        {
            return []<size_t... i>(Args && ...args, std::index_sequence<i...>)
            {
                return std::array{(i, T{std::forward<Args>(args)...})...};
            }
            (std::forward<Args>(args)..., std::make_index_sequence<n>{});
        }
    }

    template<typename... Args>
    constexpr ArrayPP(Args &&...args) : super{init(args...)}
    {}
};
