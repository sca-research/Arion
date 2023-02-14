#pragma once

#include "gadget/digest_variable_pp.hpp"
#include "gadget/field_variable.hpp"
#include "gadget/pb_variable_pp.hpp"
#include "util/array_utils.hpp"
#include <iostream>
#include <string>
#include <type_traits>

template<size_t height, typename GadHashT>
class MTreeGadget : public GadgetPP<typename GadHashT::Field>
{
public:
    using super = GadgetPP<typename GadHashT::Field>;
    using GadHash = GadHashT;
    using Field = typename GadHash::Field;

    using super::constrain;
    using super::val;

    static constexpr size_t HEIGHT = height;
    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;
    static constexpr size_t DIGEST_SIZE = GadHash::DIGEST_SIZE;
    static constexpr size_t ARITY = GadHash::BLOCK_SIZE / GadHash::DIGEST_SIZE;
    static constexpr bool HASH_ISBOOLEAN = DIGEST_SIZE < DIGEST_VARS;

    using LC = libsnark::linear_combination<Field>;
    using PbVar = PbVariablePP<Field>;
    using DigVar =
        typename std::conditional_t<HASH_ISBOOLEAN, DigestVariablePP<Field>, FieldVariable<Field>>;
    using Protoboard = libsnark::protoboard<Field>;
    using Level = std::array<DigVar, ARITY>;
    using BoolLevel = std::array<PbVar, ARITY>;

private:
    static constexpr size_t HEIGHT1 = HEIGHT - 1;

    DigVar trans;
    std::vector<Level> other;
    PbVar idx;
    std::vector<DigVar> inter;
    std::vector<Level> children;
    std::vector<GadHash> hash;
    std::vector<BoolLevel> active;

public:
    const DigVar out;

    MTreeGadget(Protoboard &pb, const DigVar &out, const DigVar &trans,
                const std::vector<Level> &other, const PbVar &idx, const std::string &ap) :
        super{pb, ap}, //
        trans{trans},  //
        other{other},  //
        idx{idx},      //
        out{out}       //
    {
        for (size_t i = 0; i < HEIGHT1; ++i)
        {
            // inputs for the hash gadget
            children.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("")));
            // active index for next level (i.e. where the previous level was output)
            active.emplace_back(make_uniform_array<BoolLevel>(pb, FMT("")));

            // result of the hash
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));

            // hash gadget
            if (i == HEIGHT1 - 1)
                hash.emplace_back(pb, children[i], out, FMT(""));
            else
                hash.emplace_back(pb, children[i], inter[i], FMT(""));
        }
    }

    void generate_r1cs_constraints()
    {
        LC sigma{0};
        Field coeff{1};

        // choices must combine to match the index
        for (size_t i = 0; i < HEIGHT1; ++i, coeff *= ARITY)
        {
            LC sigma_l{0};
            Field coeff_l{1};

            // active = 0/1, multiplied by its position k gives the level index (0/k)
            for (size_t j = 0; j < ARITY; ++j, ++coeff_l)
                sigma_l = sigma_l + active[i][j] * coeff_l;
            // multiplied by the level coefficient gives the polynomial component (k * ARITY^i)
            sigma = sigma * coeff;
        }
        constrain(sigma, 1, idx);

        // Constraints for layers
        for (size_t i = 0; i < HEIGHT1; ++i)
        {
            coeff = 0;
            // active child must equal intermediate value
            for (size_t j = 0; j < ARITY; ++j, ++coeff)
            {
                // iterate over all pieces of a single node
                for (size_t k = 0; k < DIGEST_VARS; ++k)
                {
                    // z = c ? x : y <==> z = xc + y(1 - c) <==> z - y = c(x - y)
                    if (i == 0)
                        constrain(active[i][j], trans[k] - other[i][j][k],
                                  children[i][j][k] - other[i][j][k]);
                    else
                        constrain(active[i][j], inter[i - 1][k] - other[i][j][k],
                                  children[i][j][k] - other[i][j][k]);
                }
            }
            hash[i].generate_r1cs_constraints();
        }
    }

    void generate_r1cs_witness()
    {
        // we assume idx < 2^64 (i.e. height < 64)
        // 0 ==> left, 1 ==> right
        for (size_t i = 0, uidx = val(idx).as_ulong(); i < HEIGHT1; ++i, uidx /= ARITY)
        {
            size_t rem = uidx % ARITY;

            for (size_t j = 0; j < ARITY; ++j)
            {
                if (rem == j)
                {
                    val(active[i][j]) = 1;

                    for (size_t k = 0; k < DIGEST_VARS; ++k)
                    {
                        if (i == 0)
                            val(children[i][j][k]) = val(trans[k]);
                        else
                            val(children[i][j][k]) = val(inter[i - 1][k]);
                    }
                }
                else
                {
                    val(active[i][j]) = 0;

                    for (size_t k = 0; k < DIGEST_VARS; ++k)
                        val(children[i][j][k]) = val(other[i][j][k]);
                }
            }
            hash[i].generate_r1cs_witness();
        }
    }
};
