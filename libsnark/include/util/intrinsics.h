#pragma once

#include <inttypes.h>

#ifdef _WIN32
    #include <intrin.h>
    #ifndef _bswap64
        #define _bswap64(x) _byteswap_uint64(x)
    #endif
    #ifndef _popcnt16
        #define _popcnt16(x) __popcnt16(x)
    #endif
    #ifndef _popcnt32
        #define _popcnt32(x) __popcnt(x)
    #endif
    #ifndef _popcnt64
        #define _popcnt64(x) __popcnt64(x)
    #endif
#else
    #define __int8 char
    #define __int16 short
    #define __int32 int
    #define __int64 long long

    #include <x86intrin.h>
    #define _rotr64 _lrotr
    #define _rotl64 _lrotl

    #ifdef __cplusplus
extern "C" {
    #endif
inline uint64_t _udiv128(uint64_t hi, uint64_t lo, uint64_t div, uint64_t *rem)
{
    uint64_t r, quot;

    // High bits go in RDX, low bits in RAX, quotient is in RAX, remainder is in RDX
    // for some stupid reason we must use AT&T syntax...
    asm( //
        "div %4\n"
        "mov %%rax, %0\n"
        "mov %%rdx, %1\n"
        : "=r"(quot), "=r"(r)
        : "d"(hi), "a"(lo), "r"(div));

    *rem = r;
    return quot;
}
#endif
#ifdef __cplusplus
}
#endif
