#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fq<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 1, size_t rounds = 12>
class Griffin
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_N = rounds;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;
    static constexpr size_t CIRC_N = std::min(BRANCH_N, (size_t)8);

    using Sponge = std::array<Field, BRANCH_N>;
    using CircMat = std::array<Field, CIRC_N>;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static inline CircMat circular_matrix()
    {
        if constexpr (BRANCH_N == 3)
            return CircMat{2, 1, 1};
        else if constexpr (BRANCH_N == 4)
            return CircMat{3, 2, 1, 1};
        else
        {
            static_assert(BRANCH_N % 4 == 0, "Invalid branch size");

            return CircMat{6, 4, 2, 2, 3, 2, 1, 1};
        }
    }

    static inline const Field d{5};
    static inline const Field e{inverse(d, Field{-1})};
    static inline const std::pair<Field, Field> alpha{get_irreducible_pair<Field>(~0ULL)};
    static inline const Field gamma{Field::random_element()};
    static inline const std::array<Field, ROUNDS_N * BRANCH_N> round_c{
        random_array<Field, ROUNDS_N * BRANCH_N>()};
    static inline const CircMat circ_mat{circular_matrix()};

    static void fifth(Field &x)
    {
        Field t{x};

        x *= x;
        x *= x;
        x *= t;
    }

#ifdef CURVE_ALT_BN128
    static void fifth_inv(Field &x)
    {
        // Explicit steps for 1/5 exponentiation in BN128, it is 1.4x faster than x ^= e
        Field t{x};

        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t;
        t *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, t *= t, x *= t, t *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t;
        x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, t *= t, t *= t, x *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t;
        x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, t *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t;
        t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t;
        t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, x *= t;
        t *= t, x *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t;
        t *= t, t *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, t *= t, t *= t, x *= t, t *= t;
        t *= t, x *= t, t *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t;
        t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t;
        x *= t, t *= t, t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, x *= t, t *= t;
        t *= t, x *= t, t *= t, x *= t, t *= t, t *= t, t *= t, x *= t;
    }
#else
    static void fifth_inv(Field &x)
    {
        static const auto eb = e.as_bigint();
        x ^= eb;
    }
#endif

    static void circular(Sponge &x)
    {
        if constexpr (BRANCH_N == 3)
        {
            Field old{x[0]};

            x[0] += x[1];
            x[0] += x[2];
            x[1] += x[0];
            x[2] += x[0];
            x[0] += old;
        }
        else if constexpr (BRANCH_N == 4)
        {
            Field old[2]{x[0], x[1]};

            x[0] += x[1];
            x[0] += x[2];
            x[0] += x[3];

            x[1] += x[1];
            x[1] += x[2];
            x[1] += x[0];

            x[2] += x[2];
            x[2] += x[3];
            x[2] += x[0];

            x[3] += x[3];
            x[3] += x[0];
            x[3] += old[0];
            x[0] += old[0];
            x[0] += old[0];
            x[0] += old[1];
        }
        else
        {
            static constexpr size_t BRANCH_N4 = BRANCH_N / 4;

            Sponge sum{};

            for (size_t i = 0; i < BRANCH_N4; ++i)
                for (size_t j = 0; j < BRANCH_N4; ++j)
                    for (size_t k = 0, off = 4 * (i != j); k < 4; ++k)
                        for (size_t l = 0; l < 4; ++l)
                            sum[4 * i + k] += circ_mat[off + l] * x[4 * j + l];

            x = sum;
        }
    }

    static void sbox(Sponge &x)
    {
        // Base case, y[0] = x[0]^e = x[0]^(1/d)
        fifth_inv(x[0]);

        // Base case, y[1] = x[1]^d
        fifth(x[1]);

        // Recursive case y[i] = x[i] * (L(y0,y1,old)^2 + a1*L(y0,y1,old) + a2)
        // <==> y[i] = x[i] * (L(y0,y1,old) * (L(y0,y1,old) + a1) + a2)
        // L(y1, y2, old) = gamma*y1 + y2 + old
        Field l;
        Field old{}; // old = 0 at the beginning

        for (size_t i = 2; i < BRANCH_N; ++i)
        {
            l = gamma;
            l *= x[0];
            l += x[1];
            l += old;
            old = x[i];
            x[i] = l;
            x[i] += alpha.first;
            x[i] *= l;
            x[i] += alpha.second;
            x[i] *= old;
        }
    }

    static void hash_field(Sponge &h)
    {
        // Round 0, we assume key = 0, so no key addition is ever needed
        circular(h);
        for (size_t i = 0; i < ROUNDS_N; ++i)
        {
            sbox(h);
            circular(h);
            for (size_t j = 0; j < BRANCH_N; ++j)
                h[j] += round_c[i * BRANCH_N + j];
        }
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

    Griffin() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Griffin<>;
