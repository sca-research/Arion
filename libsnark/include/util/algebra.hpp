#pragma once

#include <array>
#include <cstddef>
#include <gmpxx.h>
#include <libff/common/default_types/ec_pp.hpp>
#include <utility>

template<typename FieldT>
constexpr size_t field_size()
{
    return FieldT::num_limbs * sizeof(mp_limb_t);
}

template<typename FieldT>
mpz_class field_to_mpz(const FieldT &x)
{
    mpz_class x_mpz;

    x.as_bigint().to_mpz(x_mpz.get_mpz_t());

    return x_mpz;
}

template<typename Bigint>
mpz_class bigint_to_mpz(const Bigint &x)
{
    mpz_class x_mpz;

    x.to_mpz(x_mpz.get_mpz_t());

    return x_mpz;
}

template<typename FieldT>
int legendre(const FieldT &x)
{
    static const auto p = bigint_to_mpz(FieldT::mod);

    return mpz_legendre(field_to_mpz(x).get_mpz_t(), p.get_mpz_t());
}

template<typename FieldT>
std::pair<FieldT, FieldT> get_irreducible_pair(size_t max_trials = 32)
{
    for (size_t i = 0; i < max_trials; ++i)
    {
        std::pair<FieldT, FieldT> ab{FieldT::random_element(), FieldT::random_element()};

        if (legendre(ab.first.squared() - (ab.second + ab.second + ab.second + ab.second)) < 0)
            return ab;
    }

    return {0, 0};
}

template<typename FieldT>
FieldT inverse(const FieldT &x, const FieldT &modulus)
{
    static constexpr mp_size_t n = FieldT::num_limbs;

    using Bigint = libff::bigint<n>;

    if (x.is_zero())
        return 0;

    Bigint mod{modulus.is_zero() ? FieldT::mod : modulus.as_bigint()};
    Bigint xb{x.as_bigint()};
    Bigint v{mod}; // both source operands are destroyed by mpn_gcdext
    mp_limb_t g[n]; /* gp should have room for vn = n limbs */
    mp_limb_t s[n + 1]; /* sp should have room for vn+1 limbs */
    mp_size_t sn;

    /* computes gcd(u, v) = g = u*s + v*t, so s*u will be 1 (mod v) */
    mp_size_t gn = mpn_gcdext(g, s, &sn, xb.data, n, v.data, n);
    mp_size_t asn = std::abs(sn);

    if (!(gn == 1 && g[0] == 1))
        return 0; /* inverse exists */


    if (asn >= n)
    {
        /* if sn could require modulus reduction, do it here */
        mpn_tdiv_qr(g, xb.data, 0, s, asn, mod.data, n);
    }
    else
    {
        /* otherwise just copy it over */
        mpn_zero(xb.data, n);
        mpn_copyi(xb.data, s, asn);
    }

    /* fix up the negative sn */
    if (sn < 0 && mpn_sub_n(xb.data, mod.data, xb.data, n) != 0)
        return 0;

    return xb;
}

template<typename FieldT, size_t sz>
static std::array<FieldT, sz> random_array()
{
    std::array<FieldT, sz> a;

    std::generate(a.begin(), a.end(), []() { return FieldT::random_element(); });

    return a;
}

namespace libff
{
    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator++(Fp_model<n, m> &x)
    {
        return x += 1;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator++(Fp_model<n, m> &x, int)
    {
        auto t{x};
        ++x;

        return t;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator--(Fp_model<n, m> &x)
    {
        return x -= 1;
    }

    template<mp_size_t n, const bigint<n> &m>
    Fp_model<n, m> operator--(Fp_model<n, m> &x, int)
    {
        auto t{x};
        --x;

        return t;
    }
}; // namespace libff
