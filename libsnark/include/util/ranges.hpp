#pragma once

#if __cplusplus >= 202002L
    #include <ranges>
#else
    #include <experimental/type_traits>

template<typename T>
using IsBeginnable = decltype(std::begin(std::declval<T>()));
template<typename T>
using IsEndable = decltype(std::begin(std::declval<T>()));


template<typename T>
struct IsRange
{
    static inline constexpr bool value = std::experimental::is_detected_v<IsBeginnable, T> &&
                                         std::experimental::is_detected_v<IsEndable, T>;
};

template<typename T>
static inline constexpr bool IsRange_v = IsRange<T>::value;
#endif
