#pragma once

#ifdef _WIN32
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif
#include <cinttypes>

class Sha256
{
public:
    static constexpr size_t BLOCK_SIZE = 64;
    static constexpr size_t DIGEST_SIZE = 32;

    Sha256() = delete;

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        static constexpr uint32_t k[BLOCK_SIZE] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
            0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
            0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
            0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
            0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
            0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
            0xc67178f2,
        };

        uint32_t w[64];
        uint32_t wv[DIGEST_SIZE / sizeof(uint32_t)] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372,
                                                       0xa54ff53a, 0x510e527f, 0x9b05688c,
                                                       0x1f83d9ab, 0x5be0cd19};

        for (uint32_t i = 0; i < 16; ++i)
            w[i] = _bswap(((const uint32_t *)message)[i]);

        for (uint32_t i = 16; i < 64; ++i)
            w[i] = (_rotr(w[i - 2], 17) ^ _rotr(w[i - 2], 19) ^ w[i - 2] >> 10) + w[i - 7] +
                   (_rotr(w[i - 15], 7) ^ _rotr(w[i - 15], 18) ^ w[i - 15] >> 3) + w[i - 16];

        for (uint32_t i = 0; i < 64; ++i)
        {
            uint32_t t1 = wv[7] + (_rotr(wv[4], 6) ^ _rotr(wv[4], 11) ^ _rotr(wv[4], 25)) +
                          ((wv[4] & wv[5]) ^ (~wv[4] & wv[6])) + k[i] + w[i];

            uint32_t t2 = (_rotr(wv[0], 2) ^ _rotr(wv[0], 13) ^ _rotr(wv[0], 22)) +
                          ((wv[0] & wv[1]) ^ (wv[0] & wv[2]) ^ (wv[1] & wv[2]));

            wv[7] = wv[6];
            wv[6] = wv[5];
            wv[5] = wv[4];
            wv[4] = wv[3] + t1;
            wv[3] = wv[2];
            wv[2] = wv[1];
            wv[1] = wv[0];
            wv[0] = t1 + t2;
        }

        wv[0] += 0x6a09e667;
        wv[1] += 0xbb67ae85;
        wv[2] += 0x3c6ef372;
        wv[3] += 0xa54ff53a;
        wv[4] += 0x510e527f;
        wv[5] += 0x9b05688c;
        wv[6] += 0x1f83d9ab;
        wv[7] += 0x5be0cd19;

        for (uint32_t i = 0; i < 8; i++)
            ((uint32_t *)digest)[i] = _bswap(wv[i]);
    }

    static void hash_add(void *x, const void *y)
    {
        uint8_t *xb = (uint8_t *)x;
        const uint8_t *yb = (const uint8_t *)y;

        for (size_t i = 0; i < DIGEST_SIZE; ++i)
            xb[i] ^= yb[i];
    }
}; // namespace sha256
