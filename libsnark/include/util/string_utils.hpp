#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <ostream>
#include <random>
#include <string>


#ifdef USE_LIBFF
    #include <gmpxx.h>
    #include <libff/common/default_types/ec_pp.hpp>
#endif

#include "util/bit_pack.hpp"
#include "util/ranges.hpp"

// Convert an ASCII character to an hex digit
constexpr inline uint8_t ascii_to_digit(char c)
{
    c -= '0';
    if (c > 9)
        c -= 'A' - '0' - 10;
    if (c > 15)
        c -= 'a' - 'A';

    return c;
}

namespace detail
{
    template<size_t sz>
    struct StrToHex
    {
        std::array<uint8_t, sz / 2> v;

        constexpr StrToHex(const char (&str)[sz])
        {
            for (size_t i = 0; i < v.size(); ++i)
                v[i] = ascii_to_digit(str[i * 2]) << 4 | ascii_to_digit(str[i * 2 + 1]);
        }
    };

    template<typename T, bool rev = false>
    constexpr T str_to_t(const char *s, size_t sz)
    {
        T x = 0;

        for (size_t i = 0; i < sizeof(x) && i < sz; ++i)
            if constexpr (rev)
                x |= static_cast<decltype(x)>(s[i]) << (i * 8);
            else
                x |= static_cast<decltype(x)>(s[i]) << ((sizeof(x) - i - 1) * 8);

        return x;
    }
} // namespace detail

// Convert (at compile time) a string literal of hexadecimal digits to an std::array of
// hexadecimal values.

#define CAT0(x, y) x##y
#define CAT(x, y) CAT0(x, y)

#define STR_UTILS_STR0(x) #x
#define STR_UTILS_STR(x) STR_UTILS_STR0(x)

#define OPERATOR_STR_TO_INTTYPE_NAMED(sign, size, type, name)                                      \
    constexpr type name(const char *s, size_t sz)                                                  \
    {                                                                                              \
        return detail::str_to_t<type, false>(s, sz);                                               \
    }                                                                                              \
                                                                                                   \
    constexpr type CAT(name, r)(const char *s, size_t sz)                                          \
    {                                                                                              \
        return detail::str_to_t<type, true>(s, sz);                                                \
    }

#define OPERATOR_STR_TO_INTTYPE(sign, size)                                                        \
    OPERATOR_STR_TO_INTTYPE_NAMED(sign, size, CAT(CAT(CAT(sign, int), size), _t),                  \
                                  CAT(CAT(operator""_, sign), size))

OPERATOR_STR_TO_INTTYPE(, 64)  // ""_64, ""_64r
OPERATOR_STR_TO_INTTYPE(u, 64) // ""_u64, ""_u64r
OPERATOR_STR_TO_INTTYPE(, 32)  // ""_32, ""_32r
OPERATOR_STR_TO_INTTYPE(u, 32) // ""_u32, ""_u32r
OPERATOR_STR_TO_INTTYPE(, 16)  // ""_16, ""_16r
OPERATOR_STR_TO_INTTYPE(u, 16) // ""_u16, ""_u16r
OPERATOR_STR_TO_INTTYPE(, 8)   // ""_8, ""_8r
OPERATOR_STR_TO_INTTYPE(u, 8)  // ""_u8, ""_u8r

#if __cplusplus >= 202002L
#define BIGHEX(X) CAT(STR_UTILS_STR(X), _x)
template<detail::StrToHex cvt>
constexpr auto operator""_x()
{
    return cvt.v;
}
#else
#define BIGHEX(X) CAT(0x, CAT(X, _x))
template<char... str>
constexpr auto operator""_x()
{
    std::array<uint8_t, sizeof...(str)> s{str...};
    std::array<uint8_t, (sizeof...(str) - 2) / 2> v;

    for (size_t i = 0; i < v.size(); ++i)
        v[i] = ascii_to_digit(s[(i + 1) * 2]) << 4 | ascii_to_digit(s[(i + 1) * 2 + 1]);

    return v;
}
#endif


// Print pairs
template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &pair)
{
    return os << '(' << pair.first << ", " << pair.second << ')';
}

// Print char16_t
inline std::ostream &operator<<(std::ostream &os, const char16_t *str)
{
    while (*str)
        os << (char)*str++;

    return os;
}

// Print hexadeciaml numbers in a compact way without leaking
#define hexout(x) std::hex << std::setw(sizeof(x) * 2) << std::setfill('0') << (x) << std::dec

// Print hexadeciaml ranges in a compact way without leaking
#define hexoutit(x) std::hex << std::setw(sizeof(*std::begin(x)) * 2) << (x) << std::dec

template<typename T>
std::string hexdump(T *arr, size_t sz, bool up = false, bool rev = false, size_t space = 0)
{
    const char *TO_CHAR = up ? "0123456789ABCDEF" : "0123456789abcdef";

    space /= 2;
    size_t base_len = sz * sizeof(T) * 2;
    size_t num_spaces = (space ? base_len / space : 0) * 2;
    size_t dump_len = base_len + num_spaces;
    uint8_t *data = (uint8_t *)arr;


    std::string s(dump_len, 0);

    char *sp = s.data();

    if (rev)
        for (size_t i = 0; i < sz; ++i)
        {
            size_t k = i * sizeof(T);

            for (size_t j = 0; j < sizeof(T); ++j)
            {
                if (space && j && j % space == 0)
                    *(sp++) = ' ';

                *(sp++) = TO_CHAR[data[k + sizeof(T) - j - 1] >> 4];
                *(sp++) = TO_CHAR[data[k + sizeof(T) - j - 1] & 0xF];
            }
        }
    else
        for (size_t i = 0; i < sz * sizeof(T); ++i)
        {
            if (space && i && i % space == 0)
                *(sp++) = ' ';
            *(sp++) = TO_CHAR[data[i] >> 4];
            *(sp++) = TO_CHAR[data[i] & 0xF];
        }

    return s;
}


template<typename Iter>
std::string hexdump(const Iter &first, const Iter &last, bool up = false, bool rev = false,
                    size_t space = 0)
{
    return hexdump(&*first, std::distance(first, last), up, rev, space);
}

std::string hexdump(const std::vector<bool> &v, bool up = false, bool rev = false, size_t space = 0)
{
    std::vector<uint8_t> w(v.size() / CHAR_BIT + (v.size() % CHAR_BIT != 0));

    pack_bits(w.data(), v);

    return hexdump(w.begin(), w.end(), up, rev, space);
}

#if __cplusplus >= 202002L
template<std::ranges::range Range>
std::string hexdump(const Range &it, bool up = false, bool rev = false, size_t space = 0)
{
    return hexdump(std::ranges::begin(it), std::ranges::end(it), up, rev, space);
}
#else
template<typename T, size_t sz>
std::string hexdump(const T (&arr)[sz], bool up = false, bool rev = false, size_t space = 0)
{
    return hexdump(arr, sz, up, rev, space);
}

template<typename T>
std::string hexdump(const std::vector<T> &it, bool up = false, bool rev = false, size_t space = 0)
{
    return hexdump(std::begin(it), std::end(it), up, rev, space);
}

template<typename T, size_t sz>
std::string hexdump(const std::array<T, sz> &it, bool up = false, bool rev = false,
                    size_t space = 0)
{
    return hexdump(std::begin(it), std::end(it), up, rev, space);
}
#endif

#ifdef USE_LIBFF
template<mp_size_t limbs>
std::string hexdump(const libff::bigint<limbs> &x, bool up = false, bool rev = false,
                    size_t space = 0)
{
    mp_limb_t buff[limbs]{}; // limbs are 64-bit wide
    mpz_class tmp;

    x.to_mpz(tmp.get_mpz_t());
    mpz_export(buff, NULL, 1, 1, 0, 0, tmp.get_mpz_t());

    return hexdump(buff, up, rev, space);
}

template<typename FieldT>
std::string hexdump(const FieldT &x, bool up = false, bool rev = false, size_t space = 0)
{
    return hexdump(x.as_bigint(), up, rev, space);
}
#endif

inline std::string random_string(size_t len, std::string_view alphabet)
{
    static std::mt19937_64 rng{std::random_device{}()};

    std::uniform_int_distribution<size_t> dist{0, alphabet.size() - 1};
    std::string str(len, 0);

    std::generate(str.begin(), str.end(), [&]() { return alphabet[dist(rng)]; });

    return str;
}

#if __cplusplus >= 202002L
template<typename T>
concept beginChar = requires(T x)
{
    {
        *std::begin(x)
        } -> std::convertible_to<const char>;
};

template<typename T>
concept rangeNotChar = std::ranges::range<T> && !std::convertible_to<T, const char *> &&
                       !beginChar<T>;
#endif

// Print iterables (except std::string, std::string_view and char arrays)
#if __cplusplus >= 202002L
template<rangeNotChar Range>
std::ostream &operator<<(std::ostream &os, const Range &it)
{
    using T = decltype(*std::begin(it));
    using cast_type = std::conditional_t<std::is_convertible_v<T, unsigned char>, int, T>;

    size_t fix_width = os.width();

    os << std::setw(0) << '{';
    if (std::begin(it) != std::end(it))
    {
        auto i = std::begin(it);
        auto last = std::end(it);

        for (--last; i != last; ++i)
            if (fix_width)
                os << std::setw(fix_width) << std::setfill('0') << static_cast<cast_type>(*i)
                   << ", ";
            else
                os << static_cast<cast_type>(*i) << ", ";


        os << static_cast<cast_type>(*i);
    }

    return os << '}';
}
#else

template<typename Range,
         std::enable_if_t<IsRange_v<Range> && !std::is_same_v<Range, std::string> &&
                              !std::is_same_v<Range, std::string_view> &&
                              !std::is_same_v<decltype(*std::begin(std::declval<Range>())), char>,
                          bool> = true>
std::ostream &operator<<(std::ostream &os, const Range &it)
{
    using T = decltype(*std::begin(it));
    using cast_type = std::conditional_t<std::is_convertible_v<T, unsigned char>, int, T>;

    size_t fix_width = os.width();

    os << std::setw(0) << '{';
    if (std::begin(it) != std::end(it))
    {
        auto i = std::begin(it);
        auto last = std::end(it);

        for (--last; i != last; ++i)
            if (fix_width)
                os << std::setw(fix_width) << std::setfill('0') << static_cast<cast_type>(*i)
                   << ", ";
            else
                os << static_cast<cast_type>(*i) << ", ";


        os << static_cast<cast_type>(*i);
    }

    return os << '}';
}
#endif

#if __cplusplus >= 202002L
template<std::ranges::range Range>
constexpr void hexstring_to_range(std::string_view str, Range &out)
{
    auto out_it = std::ranges::begin(out);

    for (size_t i = 0; i < str.size(); i += 2)
        *out_it++ = ascii_to_digit(str[i]) << 4 | ascii_to_digit(str[i + 1]);
}
#endif

#undef OPERATOR_STR_TO_INTTYPE
#undef OPERATOR_STR_TO_INTTYPE_NAMED
