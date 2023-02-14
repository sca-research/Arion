#pragma once

#include "util/algebra.hpp"

template<typename FieldT = libff::Fq<libff::default_ec_pp>, size_t rate = 2, size_t capacity = 1,
         size_t rounds = 9>
class Arion
{
public:
    using Field = FieldT;

    static constexpr size_t RATE = rate;
    static constexpr size_t CAPACITY = capacity;
    static constexpr size_t ROUNDS_N = rounds;
    static constexpr size_t DIGEST_SIZE = field_size<Field>();
    static constexpr size_t BLOCK_SIZE = DIGEST_SIZE * RATE;
    static constexpr size_t BRANCH_N = RATE + CAPACITY;

    using Sponge = std::array<Field, BRANCH_N>;

    static inline const struct Init
    {
        Init() { libff::default_ec_pp::init_public_params(); }
    } init;

    static inline const Field d{5};
    static inline const Field e{inverse(d, Field{-1})};

    static inline const std::pair<Field, Field> alpha{get_irreducible_pair<Field>(~0ULL)};
    static inline const Field beta1{Field::random_element()};

    static inline const std::array<Field, ROUNDS_N * BRANCH_N> round_c{
        random_array<Field, ROUNDS_N * BRANCH_N>()};

    static Sponge iota()
    {
        Sponge s;

        std::iota(s.begin(), s.end(), Field{1});

        return s;
    }

    static inline const Sponge circ_mat{iota()};


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
        static const auto eb{e.as_bigint()};

        x ^= eb;
    }
#endif


    static void circular(Sponge &x)
    {
        if constexpr (BRANCH_N == 3)
        {
            // Explicit optimized steps for 3x3 circular matrix (it's 2.5x faster)
            Sponge s;

            s[0] = x[0];
            s[0] += x[1];
            s[2] = s[0];
            s[0] += x[2];
            s[2] += s[0];
            s[2] += x[1];
            s[0] += x[2];
            s[1] = s[0];
            s[1] += x[0];
            s[1] += x[0];
            s[0] += x[1];
            s[0] += x[2];
            x = s;
        }
        else
        {
            Field sigma{x[0]};
            Field old{x[0]};

            for (size_t i = 1; i < BRANCH_N; ++i)
                sigma += x[i];

            x[0] = sigma;
            for (size_t i = 1; i < BRANCH_N; ++i)
                x[0] += circ_mat[i - 1] * x[i];

            for (size_t i = 1; i < BRANCH_N; ++i)
            {
                std::swap(x[i], old);
                x[i] *= circ_mat[BRANCH_N - 1];
                x[i] += x[i - 1];
                x[i] -= sigma;
            }
        }
    }

    static void gtds(Sponge &x)
    {
        Sponge f;
        Field t;
        Field sigma;

        // Base case, f(x) = x[n]^e = x[n]^(1/d)
        f[BRANCH_N - 1] = x[BRANCH_N - 1];
        fifth_inv(f[BRANCH_N - 1]);

        // Recursive case: f(x) = x[i]^d * g(x) + h(x)
        for (size_t i = BRANCH_N - 2; i != (size_t)~0; --i)
        {
            // f(x[i]) = x[i]^d
            f[i] = x[i];
            fifth(f[i]);
            // sigma = sum_{j=i+1}^{BRANCH_N}{x[j] + f[j]}
            sigma = x[i + 1] + f[i + 1];
            for (size_t j = i + 2; j < BRANCH_N; ++j)
                sigma += x[j] + f[j];

            // t = g(x) = sigma^2 + alpha1*sigma + alpha2
            t = sigma;
            t += alpha.first;
            t *= sigma;
            t += alpha.second;
            f[i] *= t;

            // t = h(x) = sigma^2 + beta1*sigma
            t = sigma;
            t += beta1;
            t *= sigma;
            f[i] += t;
        }

        x = f;
    }

    static void hash_field(Sponge &h)
    {
        // Round 0, we assume key = 0, so no key addition is ever needed
        circular(h);

        for (size_t i = 0; i < ROUNDS_N; ++i)
        {
            gtds(h);
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

    Arion() = delete;
};

// Valid realizations of a template class must be initialized before main()!
template class Arion<>;
