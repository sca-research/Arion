#pragma once

#include <climits>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "util/ranges.hpp"

static inline void pack_bits(uint8_t *v, const std::vector<bool> &bv)
{
    for (size_t i = 0, sz = bv.size() / CHAR_BIT; i < sz; ++i)
    {
        v[i] = 0;
        v[i] |= bv[i * 8 + 0] << 7;
        v[i] |= bv[i * 8 + 1] << 6;
        v[i] |= bv[i * 8 + 2] << 5;
        v[i] |= bv[i * 8 + 3] << 4;
        v[i] |= bv[i * 8 + 4] << 3;
        v[i] |= bv[i * 8 + 5] << 2;
        v[i] |= bv[i * 8 + 6] << 1;
        v[i] |= bv[i * 8 + 7] << 0;
    }
}

static inline void unpack_bits(std::vector<bool> &bv, const void *data)
{
    const uint8_t *v = static_cast<const uint8_t *>(data);

    for (size_t i = 0, sz = bv.size() / CHAR_BIT; i < sz; ++i)
    {
        bv[i * 8 + 0] = v[i] >> 7 & 1;
        bv[i * 8 + 1] = v[i] >> 6 & 1;
        bv[i * 8 + 2] = v[i] >> 5 & 1;
        bv[i * 8 + 3] = v[i] >> 4 & 1;
        bv[i * 8 + 4] = v[i] >> 3 & 1;
        bv[i * 8 + 5] = v[i] >> 2 & 1;
        bv[i * 8 + 6] = v[i] >> 1 & 1;
        bv[i * 8 + 7] = v[i] >> 0 & 1;
    }
}

static inline std::vector<bool> unpack_bits(const void *data, size_t sz)
{
    std::vector<bool> bv(sz * CHAR_BIT);

    unpack_bits(bv, data);

    return bv;
}

template<typename It>
static inline std::vector<bool> unpack_bits(It first, It last)
{
    return unpack_bits(&*first, std::distance(first, last));
}
#if __cplusplus >= 202002L
template<std::ranges::range Range>
static inline std::vector<bool> unpack_bits(const Range &r)
{
    return unpack_bits(std::ranges::begin(r), std::ranges::end(r));
}
#else
template<typename Range, std::enable_if_t<IsRange_v<Range>, bool> = true>
static inline std::vector<bool> unpack_bits(const Range &r)
{
    return unpack_bits(std::data(r), std::size(r));
}
#endif
