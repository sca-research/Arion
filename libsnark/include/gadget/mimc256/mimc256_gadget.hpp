#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"
#include "hash/mimc256.hpp"

template<typename FieldT>
class mimc256_two_to_one_hash_gadget : public GadgetPP<FieldT>
{
public:
    using super = GadgetPP<FieldT>;
    using Hash = Mimc256<FieldT>;
    using DigVar = FieldVariable<FieldT>;
    using LC = libsnark::linear_combination<FieldT>;
    using PbVar = libsnark::pb_variable<FieldT>;
    using Field = FieldT;

    using super::constrain;
    using super::val;

    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    static inline const auto &round_c = Hash::round_c;

    const DigVar x;
    const DigVar y;
    const DigVar out;

private:
    static constexpr size_t INTER_N = 2 + 2 * (ROUNDS_N - 1) + 2 + 2 * (ROUNDS_N - 2) + 1;

    std::vector<PbVar> inter;

public:
    mimc256_two_to_one_hash_gadget(libsnark::protoboard<FieldT> &pb, const DigVar &x,
                                   const DigVar &y, const DigVar &out,
                                   const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}
    {
        for (size_t i = 0; i < INTER_N; ++i)
        {
            inter.emplace_back();
            inter.back().allocate(pb, FMT(annotation_prefix, "_mimc256_inter_%llu", i));
        }
    }

    void generate_r1cs_constraints()
    {
        size_t i = 0;
        LC t;

        // x^2
        i += constrain(x[0], x[0], inter[i]);

        // x^3
        i += constrain(inter[i - 1], x[0], inter[i]);

        // Loop body: (x + c)^3
        for (size_t j = 0; j < ROUNDS_N - 1; ++j)
        {
            t = inter[i - 1] + round_c[j];
            i += constrain(t, t, inter[i]);
            i += constrain(inter[i - 1], t, inter[i]);
        }

        // (x + y)^3
        t = inter[i - 1] + y[0];
        i += constrain(t, t, inter[i]);
        i += constrain(inter[i - 1], t, inter[i]);

        // Loop body: x = (x + y + c)^3
        for (size_t j = 0; j < ROUNDS_N - 2; ++j)
        {
            t = inter[i - 1] + y[0] + round_c[j];
            i += constrain(t, t, inter[i]);
            i += constrain(inter[i - 1], t, inter[i]);
        }
        t = inter[i - 1] + y[0] + round_c[ROUNDS_N - 2];
        constrain(t, t, inter[i]);
        constrain(inter[i], t, out[0] - y[0]);
    }

    void generate_r1cs_witness()
    {
        size_t i = 0;
        FieldT t;

        // x^2
        val(inter[i]) = val(x[0]) * val(x[0]);
        ++i;
        // x^3
        val(inter[i]) = val(inter[i - 1]) * val(x[0]);
        ++i;

        // Loop body: (x + c)^3
        for (size_t j = 0; j < ROUNDS_N - 1; ++j)
        {
            // (x + c)
            t = val(inter[i - 1]) + round_c[j];
            val(inter[i]) = t * t;
            ++i;

            // ^3
            val(inter[i]) = val(inter[i - 1]) * t;
            ++i;
        }

        // (x + y)^2
        t = val(inter[i - 1]) + val(y[0]);
        val(inter[i]) = t * t;
        ++i;

        // ^3
        val(inter[i]) = val(inter[i - 1]) * t;
        ++i;

        // Loop body: x = (x + y + c)^3
        for (size_t j = 0; j < ROUNDS_N - 2; ++j)
        {
            t = val(inter[i - 1]) + val(y[0]) + round_c[j];
            // (x + y + c)^2
            val(inter[i]) = t * t;
            ++i;
            // ^3
            val(inter[i]) = val(inter[i - 1]) * t;
            ++i;
        }
        // (x + y + c)^2
        t = val(inter[i - 1]) + val(y[0]) + round_c[ROUNDS_N - 2];
        val(inter[i]) = t * t;
        // ^3 + y
        val(out[0]) = val(inter[i]) * t;
        val(out[0]) += val(y[0]);
    }
};
