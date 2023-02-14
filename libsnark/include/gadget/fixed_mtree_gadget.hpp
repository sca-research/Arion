#pragma once

#include "gadget/field_variable.hpp"
#include <iostream>
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>
#include <string>
#include <type_traits>

template<typename GadHash>
class FixedMTreeGadget : public libsnark::gadget<typename GadHash::Field>
{
public:
    using FieldT = typename GadHash::Field;
    using super = libsnark::gadget<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;
    static constexpr size_t DIGEST_SIZE = GadHash::DIGEST_SIZE;
    static constexpr bool HASH_ISBOOLEAN = DIGEST_SIZE < DIGEST_VARS;

    using DigVar = typename std::conditional_t<HASH_ISBOOLEAN, libsnark::digest_variable<FieldT>,
                                               FieldVariable<FieldT>>;
    using Protoboard = libsnark::protoboard<FieldT>;

private:
    DigVar trans;
    std::vector<DigVar> other;
    std::vector<DigVar> inter;
    std::vector<GadHash> foo_hash;
    size_t height;
    size_t idx;

public:
    const DigVar out;

    FixedMTreeGadget(Protoboard &pb, const DigVar &out, const DigVar &trans,
                     const std::vector<DigVar> &other, size_t trans_idx, const std::string &ap) :
        super{pb, ap},
        trans{trans}, other{other}, height{other.size() + 1}, idx{trans_idx}, out{out}
    {
        inter.emplace_back(pb, DIGEST_VARS, FMT(""));
        if (trans_idx & 1)
            foo_hash.emplace_back(pb, other[0], trans, inter[0], FMT(""));
        else
            foo_hash.emplace_back(pb, trans, other[0], inter[0], FMT(""));
        trans_idx >>= 1;

        for (size_t i = 1; i < height - 2; ++i, trans_idx >>= 1)
        {
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, other[i], inter[i - 1], inter[i], FMT(""));
            else
                foo_hash.emplace_back(pb, inter[i - 1], other[i], inter[i], FMT(""));
        }
        if (trans_idx & 1)
            foo_hash.emplace_back(pb, other[height - 2], inter[height - 3], out, FMT(""));
        else
            foo_hash.emplace_back(pb, inter[height - 3], other[height - 2], out, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        for (auto &&x : foo_hash)
            x.generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        for (auto &&x : foo_hash)
            x.generate_r1cs_witness();
    }
};
