#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fq<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 1,
         size_t rounds_f = 4, size_t rounds_p = 57>
class Poseidon5
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_f_N = rounds_f;
    static constexpr size_t ROUNDS_P_N = rounds_p;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;
    static constexpr size_t ROUNDS_F_N = 2 * ROUNDS_f_N;
    static constexpr size_t ROUNDS_N = ROUNDS_F_N + ROUNDS_P_N;
    static constexpr size_t CONST_N = BRANCH_N * ROUNDS_N;

    using Sponge = std::array<Field, BRANCH_N>;


    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    // Not an actual MDS, only use for benchmarks!
    static inline std::array<Field, BRANCH_N * BRANCH_N> mds_matrix()
    {
        return random_array<Field, BRANCH_N * BRANCH_N>();
    }

    static inline const std::array<Field, CONST_N> round_c{random_array<Field, CONST_N>()};
    static inline const std::array<Field, BRANCH_N * BRANCH_N> mds_mat{mds_matrix()};

    static void fifth(Field &x)
    {
        Field t{x};

        x *= x;
        x *= x;
        x *= t;
    }

    static void matmul(Sponge &arr)
    {
        Sponge sum;

        for (size_t i = 0; i < BRANCH_N; ++i)
        {
            sum[i] = mds_mat[i * BRANCH_N] * arr[0];
            for (size_t j = 1; j < BRANCH_N; ++j)
                sum[i] += mds_mat[i * BRANCH_N + j] * arr[j];
        }

        arr = sum;
    }


    static Field hash_field(Sponge &h)
    {
        size_t c = 0;

        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[c++];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(h[j]);

            matmul(h);
        }

        for (size_t i = 0; i < ROUNDS_P_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[c++];

            fifth(h[0]);
            matmul(h);
        }

        for (size_t i = 0; i < ROUNDS_f_N; ++i)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[c++];

            for (size_t j = 0; j < BRANCH_N; ++j)
                fifth(h[j]);

            matmul(h);
        }

        return h[0];
    }

    static void hash_oneblock(uint8_t *digest, const void *message)
    {
        Sponge h{};
        mpz_class tmp;

        for (size_t i = 0; i < RATE; ++i)
        {
            mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0,
                       (const char *)message + DIGEST_SIZE * i);

            h[i] = Field{tmp.get_mpz_t()};
        }

        hash_field(h);
        h[0].as_bigint().to_mpz(tmp.get_mpz_t());

        memset(digest, 0, DIGEST_SIZE);
        mpz_export(digest, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }

    static void hash_add(void *x, const void *y)
    {
        mpz_class tmp;

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, x);
        Field xf{tmp.get_mpz_t()};

        mpz_import(tmp.get_mpz_t(), DIGEST_SIZE, 1, 1, 0, 0, y);
        Field yf{tmp.get_mpz_t()};

        xf += yf;

        xf.as_bigint().to_mpz(tmp.get_mpz_t());
        memset(x, 0, DIGEST_SIZE);
        mpz_export(x, NULL, 1, 1, 0, 0, tmp.get_mpz_t());
    }

    Poseidon5() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Poseidon5<>;