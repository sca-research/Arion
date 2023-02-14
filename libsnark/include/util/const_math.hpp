#pragma once

template<typename T>
constexpr T pow(T x, T e)
{
    // Square and multiply algorithm
    T r = 1;

    while (e)
    {
        if (e & 1)
            r *= x;
        x *= x;
        e >>= 1;
    }

    return r;
}

template<typename T>
constexpr T pow_sum(T x, T s, T e)
{
    // Computes sum_{i = s}^{e-1} x^i
    return (pow(x, e) - pow(x, s)) / (x - 1);
}
