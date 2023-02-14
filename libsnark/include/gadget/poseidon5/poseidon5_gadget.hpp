#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"
#include "gadget/pb_variable_pp.hpp"
#include "hash/poseidon5.hpp"

template<typename Poseidon5>
class Poseidon5Gadget : public GadgetPP<typename Poseidon5::Field>
{
public:
    using Field = typename Poseidon5::Field;
    using Hash = Poseidon5;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;

    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ROUNDS_f = Hash::ROUNDS_f_N;
    static constexpr size_t ROUNDS_F = Hash::ROUNDS_F_N;
    static constexpr size_t ROUNDS_P = Hash::ROUNDS_P_N;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename Poseidon5::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER0_N = 3 * ROUNDS_N;
    static constexpr size_t INTERk_N = 3 * ROUNDS_F;

    static constexpr auto &rc = Hash::round_c;
    static constexpr auto &mds = Hash::mds_mat;

    std::vector<PbVariablePP<Field>> inter[BRANCH_N];

public:
    Poseidon5Gadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                    const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        in{in}, out{out}
    {
        for (size_t i = 0; i < INTER0_N; ++i)
            inter[0].emplace_back(pb, FMT(""));

        for (size_t i = 1; i < BRANCH_N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        LC s[BRANCH_N];
        LC t[BRANCH_N]{};
        size_t i[BRANCH_N]{};
        size_t ri = 0;

        for (size_t j = 0; j < RATE; ++j)
            t[j] = in[j][0];

        // INITIAL FULL SBOX
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // (x+c)^5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = t[k] + rc[ri++];
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }
            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = mds[k * BRANCH_N] * inter[0][i[0] - 1];
                for (size_t l = 1; l < BRANCH_N; ++l)
                    t[k] = t[k] + mds[k * BRANCH_N + l] * inter[l][i[l] - 1];
            }
        }
        // PARTIAL SBOX
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // (x+c)^5 for first element
            t[0] = t[0] + rc[ri++];
            i[0] += constrain(t[0], t[0], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], inter[0][i[0] - 1], inter[0][i[0]]);
            i[0] += constrain(inter[0][i[0] - 1], t[0], inter[0][i[0]]);

            // MDS multiplication, other elements get only addition
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = mds[k * BRANCH_N] * inter[0][i[0] - 1];
                for (size_t l = 1; l < BRANCH_N; ++l)
                    s[k] = s[k] + mds[k * BRANCH_N + l] * (t[l] + rc[ri + l - 1]);
            }
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];

            ri += BRANCH_N - 1;
        }

        // FINAL FULL SBOX
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // (x+c)^5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = t[k] + rc[ri++];
                i[k] += constrain(t[k], t[k], inter[k][i[k]]);
                i[k] += constrain(inter[k][i[k] - 1], inter[k][i[k] - 1], inter[k][i[k]]);
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }
            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = mds[k * BRANCH_N] * inter[0][i[0] - 1];
                for (size_t l = 1; l < BRANCH_N; ++l)
                    t[k] = t[k] + mds[k * BRANCH_N + l] * inter[l][i[l] - 1];
            }
        }
        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        Field s[BRANCH_N];
        Field t[BRANCH_N]{};
        size_t i[BRANCH_N]{};
        size_t ri = 0;

        for (size_t j = 0; j < RATE; ++j)
            t[j] = val(in[j][0]);

        // INITIAL FULL SBOX
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // (x+c)^5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] += rc[ri++];
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = t[k] * val(inter[k][i[k] - 1]);
                ++i[k];
            }
            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = mds[k * BRANCH_N] * val(inter[0][i[0] - 1]);
                for (size_t l = 1; l < BRANCH_N; ++l)
                    t[k] += mds[k * BRANCH_N + l] * val(inter[l][i[l] - 1]);
            }
        }

        // PARTIAL SBOX
        for (size_t j = 0; j < ROUNDS_P; ++j)
        {
            // (x+c)^5 only for first element
            t[0] += rc[ri++];
            val(inter[0][i[0]]) = t[0] * t[0];
            ++i[0];
            val(inter[0][i[0]]) = val(inter[0][i[0] - 1]) * val(inter[0][i[0] - 1]);
            ++i[0];
            val(inter[0][i[0]]) = t[0] * val(inter[0][i[0] - 1]);
            ++i[0];

            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                s[k] = mds[k * BRANCH_N] * val(inter[0][i[0] - 1]);
                for (size_t l = 1; l < BRANCH_N; ++l)
                    s[k] += mds[k * BRANCH_N + l] * (t[l] + rc[ri + l - 1]);
            }
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = s[k];

            ri += BRANCH_N - 1;
        }

        // FINAL FULL SBOX
        for (size_t j = 0; j < ROUNDS_f; ++j)
        {
            // (x+c)^5
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] += rc[ri++];
                val(inter[k][i[k]]) = t[k] * t[k];
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]) * val(inter[k][i[k] - 1]);
                ++i[k];
                val(inter[k][i[k]]) = t[k] * val(inter[k][i[k] - 1]);
                ++i[k];
            }
            // MDS multiplication
            for (size_t k = 0; k < BRANCH_N; ++k)
            {
                t[k] = mds[k * BRANCH_N] * val(inter[0][i[0] - 1]);
                for (size_t l = 1; l < BRANCH_N; ++l)
                    t[k] += mds[k * BRANCH_N + l] * val(inter[l][i[l] - 1]);
            }
        }

        val(out[0]) = t[0];
    }
};
