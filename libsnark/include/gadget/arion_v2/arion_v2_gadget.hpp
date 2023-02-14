#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"
#include "gadget/pb_variable_pp.hpp"
#include "hash/arion_v2.hpp"

template<typename ArionV2>
class ArionV2Gadget : public GadgetPP<typename ArionV2::Field>
{
public:
    using Field = typename ArionV2::Field;
    using Hash = ArionV2;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;


    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename ArionV2::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t N = BRANCH_N - 1;
    static constexpr size_t INTERn_N = 9 * ROUNDS_N;
    static constexpr size_t INTERk_N = 5 * ROUNDS_N;

    std::array<std::vector<PbVariablePP<Field>>, BRANCH_N> inter;

public:
    ArionV2Gadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                 const std::string &ap) :
        super{pb, ap},
        in{in}, out{out}, inter{}
    {
        for (size_t i = 0; i < INTERn_N; ++i)
            inter[N].emplace_back(pb, FMT(""));

        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        std::array<size_t, BRANCH_N> i{};
        std::array<LC, BRANCH_N> t{};
        LC sigma;
        LC g;

        for (size_t i = 0; i < RATE; ++i)
            t[i] = sigma + in[i][0];

        // FIRST CIRCULANT MATRIX
        sigma = t[0];
        g = t[0];
        for (size_t i = 1; i < BRANCH_N; ++i)
            sigma = sigma + t[i];

        t[0] = sigma;
        for (size_t i = 1; i < BRANCH_N; ++i)
            t[0] = t[0] + Hash::circ_mat[i - 1] * t[i];

        for (size_t i = 1; i < BRANCH_N; ++i)
        {
            std::swap(t[i], g);
            t[i] = Hash::circ_mat[BRANCH_N - 1] * t[i];
            t[i] = t[i] + t[i - 1];
            t[i] = t[i] - sigma;
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // y = x^(1/257) <==> y^257 = x
            // y^256 * y = x
            i[N] += constrain(inter[N][i[N]], inter[N][i[N] + 8], t[N]);
            // y^128 * y^128 = y^256
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^64 * y^64 = y^128
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^32 * y^32 = y^64
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^16 * y^16 = y^32
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^8 * y^8 = y^16
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^4 * y^4 = y^8
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^2 * y^2 = y^4
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);
            // y^1 * y^1 = y^2
            i[N] += constrain(inter[N][i[N]], inter[N][i[N]], inter[N][i[N] - 1]);

            // x^d * g(x) + h(x)
            for (size_t k = BRANCH_N - 2; k != (size_t)~0; --k)
            {
                // x^5
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], t[k], inter[k][i[k]]);

                // sigma = sum_{l=k+1}^{BRANCH_N}{x[l] + f[l]}
                sigma = t[k + 1] + inter[k + 1][i[k + 1] - 1];
                for (size_t l = k + 2; l < BRANCH_N; ++l)
                    sigma = sigma + t[l] + inter[l][i[l] - 1];

                i[k] += constrain(sigma, sigma, inter[k][i[k]]);

                // g(x) = s^2 + a1*s + a2
                g = inter[k][i[k] - 1] + Hash::alpha.first * sigma + Hash::alpha.second;
                // h(x) = s^2 + b*s
                sigma = inter[k][i[k] - 1] + Hash::beta1 * sigma;
                // y = x^d * g(x) + h(x) <==> y - h(x) = x^d * g(x)
                i[k] += constrain(inter[k][i[k] - 2], g, inter[k][i[k]] - sigma);
            }

            // CIRCULANT MATRIX
            sigma = inter[0][i[0] - 1];
            for (size_t k = 1; k < BRANCH_N; ++k)
                sigma = sigma + inter[k][i[k] - 1];

            t[0] = sigma;
            for (size_t k = 1; k < BRANCH_N; ++k)
                t[0] = t[0] + Hash::circ_mat[k - 1] * inter[k][i[k] - 1];

            for (size_t k = 1; k < BRANCH_N; ++k)
            {
                t[k] = Hash::circ_mat[BRANCH_N - 1] * inter[k - 1][i[k - 1] - 1];
                t[k] = t[k] + t[k - 1];
                t[k] = t[k] - sigma;
            }

            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + Hash::round_c[j * BRANCH_N + k];
        }

        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        std::array<size_t, BRANCH_N> i{};
        std::array<Field, BRANCH_N> t{};
        Field sigma;
        Field g;

        for (size_t i = 0; i < RATE; ++i)
            t[i] = val(in[i][0]);

        // FIRST CIRCULANT MATRIX
        sigma = t[0];
        g = t[0];
        for (size_t i = 1; i < BRANCH_N; ++i)
            sigma += t[i];

        t[0] = sigma;
        for (size_t i = 1; i < BRANCH_N; ++i)
            t[0] += Hash::circ_mat[i - 1] * t[i];

        for (size_t i = 1; i < BRANCH_N; ++i)
        {
            std::swap(t[i], g);
            t[i] *= Hash::circ_mat[BRANCH_N - 1];
            t[i] += t[i - 1];
            t[i] -= sigma;
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // y = x^(1/257)
            val(inter[N][i[N] + 8]) = t[N];
            Hash::pow_e(val(inter[N][i[N] + 8]));
            // y^2
            val(inter[N][i[N] + 7]) = val(inter[N][i[N] + 8]) * val(inter[N][i[N] + 8]);
            // y^4
            val(inter[N][i[N] + 6]) = val(inter[N][i[N] + 7]) * val(inter[N][i[N] + 7]);
            // y^8
            val(inter[N][i[N] + 5]) = val(inter[N][i[N] + 6]) * val(inter[N][i[N] + 6]);
            // y^16
            val(inter[N][i[N] + 4]) = val(inter[N][i[N] + 5]) * val(inter[N][i[N] + 5]);
            // y^32
            val(inter[N][i[N] + 3]) = val(inter[N][i[N] + 4]) * val(inter[N][i[N] + 4]);
            // y^64
            val(inter[N][i[N] + 2]) = val(inter[N][i[N] + 3]) * val(inter[N][i[N] + 3]);
            // y^128
            val(inter[N][i[N] + 1]) = val(inter[N][i[N] + 2]) * val(inter[N][i[N] + 2]);
            // y^256
            val(inter[N][i[N]]) = val(inter[N][i[N] + 1]) * val(inter[N][i[N] + 1]);
            i[N] += 9;

            // x^d * g(x) + h(x)
            for (size_t k = BRANCH_N - 2; k != (size_t)~0; --k)
            {
                // x^5
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * t[k];
                ++i[k];

                // sigma = sum_{l=k+1}^{BRANCH_N}{x[l] + f[l]}
                sigma = t[k + 1] + val(inter[k + 1][i[k + 1] - 1]);
                for (size_t l = k + 2; l < BRANCH_N; ++l)
                    sigma += t[l] + val(inter[l][i[l] - 1]);

                val(inter[k][i[k]]) = sigma * sigma;
                ++i[k];
                // g(x) = s^2 + a1*s + a2
                g = val(inter[k][i[k] - 1]) + Hash::alpha.first * sigma + Hash::alpha.second;
                // h(x) = s^2 + b*s
                sigma = val(inter[k][i[k] - 1]) + Hash::beta1 * sigma;
                // y = x^d * g(x) + h(x)
                val(inter[k][i[k]]) = val(inter[k][i[k] - 2]) * g + sigma;
                ++i[k];
            }
            // CIRCULANT MATRIX
            sigma = val(inter[0][i[0] - 1]);
            for (size_t k = 1; k < BRANCH_N; ++k)
                sigma += val(inter[k][i[k] - 1]);

            t[0] = sigma;
            for (size_t k = 1; k < BRANCH_N; ++k)
                t[0] += Hash::circ_mat[k - 1] * val(inter[k][i[k] - 1]);

            for (size_t k = 1; k < BRANCH_N; ++k)
            {
                t[k] = Hash::circ_mat[BRANCH_N - 1] * val(inter[k - 1][i[k - 1] - 1]);
                t[k] += t[k - 1];
                t[k] -= sigma;
            }

            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += Hash::round_c[j * BRANCH_N + k];
        }

        val(out[0]) = t[0];
    }

    size_t get_constraints_size() const { return N * INTERk_N + INTERn_N; }
};

template<typename FieldT>
class arionv2_two_to_one_hash_gadget : public ArionV2Gadget<ArionV2<FieldT, 2, 1>>
{
public:
    using super = ArionV2Gadget<ArionV2<FieldT, 2, 1>>;
    using Hash = typename super::Hash;
    using DigVar = typename super::DigVar;
    using Field = typename super::Field;
    using BlockVar = typename super::BlockVar;

    arionv2_two_to_one_hash_gadget(libsnark::protoboard<Field> &pb, const DigVar &x, const DigVar &y,
                                  const DigVar &out, const std::string &ap) :
        super{pb, BlockVar{x, y}, out, ap}
    {}
};
