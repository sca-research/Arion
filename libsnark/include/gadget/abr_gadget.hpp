#pragma once

#include "gadget/add_gadget.hpp"
#include "gadget/field_variable.hpp"
#include "gadget/xor_gadget.hpp"
#include "util/bit_pack.hpp"
#include "util/string_utils.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>
#include <libsnark/gadgetlib1/protoboard.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include <iostream>
#include <string>
#include <type_traits>

template<typename FieldT, typename GadHash>
class ABR_Gadget : public libsnark::gadget<FieldT>
{
public:
    using super = libsnark::gadget<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;
    static constexpr size_t DIGEST_SIZE = GadHash::DIGEST_SIZE;
    static constexpr bool HASH_ISBOOLEAN = DIGEST_SIZE < DIGEST_VARS;

    using Protoboard = libsnark::protoboard<FieldT>;
    using GadXor =
        typename std::conditional_t<HASH_ISBOOLEAN, LongXORGadget<FieldT>, LongAddGadget<FieldT>>;
    using DigVar = typename std::conditional_t<HASH_ISBOOLEAN, libsnark::digest_variable<FieldT>,
                                               FieldVariable<FieldT>>;

private:
    DigVar trans;
    DigVar other;
    std::vector<DigVar> boot;
    std::vector<DigVar> middle;
    std::vector<DigVar> otherx;
    std::vector<DigVar> otherxm;
    std::vector<DigVar> inter;
    std::vector<DigVar> interx;
    std::vector<DigVar> interxm;
    std::vector<GadHash> foo_hash;
    std::vector<GadXor> foo_xor;
    size_t height;
    size_t sub_height;
    size_t leaves_n;
    size_t trans_idx;
    bool is_leaf;

public:
    const DigVar out;

    static void xor_digest(uint8_t *a, const uint8_t *b)
    {
        for (size_t i = 0; i < GadHash::DIGEST_SIZE; ++i)
            a[i] ^= b[i];
    }

    ABR_Gadget(Protoboard &pb,                                                       //
               const DigVar &out, const DigVar &trans, const DigVar &other,          //
               const std::vector<DigVar> &middle, const std::vector<DigVar> &otherx, //
               size_t trans_idx, size_t height,                                      //
               const std::string &ap) :
        super{pb, ap},                  //
        trans{trans},                   //
        other{other},                   //
        middle{middle},                 //
        otherx{otherx},                 //
        height{height},                 //
        sub_height{height},             //
        leaves_n{1ULL << (height - 1)}, //
        trans_idx{trans_idx},           //
        is_leaf{trans_idx < leaves_n},  //
        out{out}
    {
        // if leaf, the first hash has no middle and no remainer, and the second has no remainder
        if (is_leaf)
        {
            // first hash without xor
            boot.emplace_back(pb, DIGEST_VARS, FMT(""));
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, other, trans, boot[0], FMT(""));
            else
                foo_hash.emplace_back(pb, trans, other, boot[0], FMT(""));

            trans_idx >>= 1;

            // xor middle node with left/right (xor is commutative)
            interxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            foo_xor.emplace_back(pb, boot[0], middle[0], interxm[0], FMT(""));
            otherxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            foo_xor.emplace_back(pb, middle[0], otherx[0], otherxm[0], FMT(""));

            // hash result and add remainder
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));
            interx.emplace_back(pb, DIGEST_VARS, FMT(""));
            if (trans_idx & 1)
            {
                foo_hash.emplace_back(pb, otherxm[0], interxm[0], inter[0], FMT(""));
                // the node was a right leaf, the remainder is boot[0]
                foo_xor.emplace_back(pb, inter[0], boot[0], interx[0], FMT(""));
            }
            else
            {
                foo_hash.emplace_back(pb, interxm[0], otherxm[0], inter[0], FMT(""));
                // the node was a left leaf or a middle, the remainder is otherx[0]
                foo_xor.emplace_back(pb, inter[0], otherx[0], interx[0], FMT(""));
            }
        }
        // if middle, the first hash has no remainder
        else
        {
            // translate trans_idx to be used as a path along the tree
            trans_idx >>= 1;
            do
            {
                trans_idx -= leaves_n;
                leaves_n /= 2;
                --sub_height; // it's like the tree starts at the provided middle node
            } while (trans_idx >= leaves_n);

            // by convention, other is the left node, and otherx[0] is the right node
            // xor middle node with left/right (xor is commutative)
            interxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            foo_xor.emplace_back(pb, other, trans, interxm[0], FMT(""));
            otherxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            foo_xor.emplace_back(pb, trans, otherx[0], otherxm[0], FMT(""));

            // hash result
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));
            interx.emplace_back(pb, DIGEST_VARS, FMT(""));
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, otherxm[0], interxm[0], inter[0], FMT(""));
            else
                foo_hash.emplace_back(pb, interxm[0], otherxm[0], inter[0], FMT(""));

            // add remainder
            foo_xor.emplace_back(pb, inter[0], otherx[0], interx[0], FMT(""));
        }
        trans_idx >>= 1;

        /* Given height=h, we have:
        1. h levels in the tree
        2. h - 1 hash outputs (i.e. non-leaves)
        3. h - 2 intermediate hash output (i.e. non-root)
        4. If our node is a leaf, we already computed the first two hashes (boot and interx[0])
        5. If our node is internal, we already computed the first hash (interx[0])
        Hence, we need to insert rem = (h - 2 - is_leaf) blocks.
        E.g., if our node is a leaf and height = 4, this loop won't run at all
        */
        for (size_t i = 1, rem = sub_height - 2 - is_leaf; i < rem; ++i, trans_idx >>= 1)
        {
            interxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            otherxm.emplace_back(pb, DIGEST_VARS, FMT(""));
            inter.emplace_back(pb, DIGEST_VARS, FMT(""));
            interx.emplace_back(pb, DIGEST_VARS, FMT(""));

            // xor middle node with left/right (xor is commutative)
            foo_xor.emplace_back(pb, interx[i - 1], middle[i], interxm[i], FMT(""));
            foo_xor.emplace_back(pb, otherx[i], middle[i], otherxm[i], FMT(""));

            // hash and add remainder
            if (trans_idx & 1)
            {
                foo_hash.emplace_back(pb, otherxm[i], interxm[i], inter[i], FMT(""));
                // the node was a right leaf, the remainder is interx[i-1]
                foo_xor.emplace_back(pb, inter[i], interx[i - 1], interx[i], FMT(""));
            }
            else
            {
                foo_hash.emplace_back(pb, interxm[i], otherxm[i], inter[i], FMT(""));
                // the node was a left leaf or a middle, the remainder is otherx[i]
                foo_xor.emplace_back(pb, inter[i], otherx[i], interx[i], FMT(""));
            }
        }

        // last layer must put result into "out"
        interxm.emplace_back(pb, DIGEST_VARS, FMT(""));
        otherxm.emplace_back(pb, DIGEST_VARS, FMT(""));
        inter.emplace_back(pb, DIGEST_VARS, FMT(""));

        // xor middle node with left/right (xor is commutative)
        foo_xor.emplace_back(pb, interx.back(), middle.back(), interxm.back(), FMT(""));
        foo_xor.emplace_back(pb, otherx.back(), middle.back(), otherxm.back(), FMT(""));

        // hash and add remainder
        if (trans_idx & 1)
        {
            foo_hash.emplace_back(pb, otherxm.back(), interxm.back(), inter.back(), FMT(""));
            // the node was a right leaf, the remainder is interx[i-1]
            foo_xor.emplace_back(pb, inter.back(), interx.back(), out, FMT(""));
        }
        else
        {
            foo_hash.emplace_back(pb, interxm.back(), otherxm.back(), inter.back(), FMT(""));
            // the node was a left leaf or a middle, the remainder is otherx[i]
            foo_xor.emplace_back(pb, inter.back(), otherx.back(), out, FMT(""));
        }
    }

    void generate_r1cs_constraints()
    {
        for (auto &&x : foo_xor)
            x.generate_r1cs_constraints();

        for (auto &&x : foo_hash)
            x.generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        size_t i_h = 0;
        size_t i_x = 0;

        // generate the "boot" node
        if (is_leaf)
            foo_hash[i_h++].generate_r1cs_witness();

        // generate the "interx" nodes and the output node
        for (size_t i = 0; i < sub_height - 2; ++i)
        {
            foo_xor[i_x++].generate_r1cs_witness();
            foo_xor[i_x++].generate_r1cs_witness();
            foo_hash[i_h++].generate_r1cs_witness();
            foo_xor[i_x++].generate_r1cs_witness();
        }
    }
};
