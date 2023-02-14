#pragma once

#include "util/bit_pack.hpp"
#include "xor_gadget.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>
#include <libsnark/gadgetlib1/protoboard.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include <iostream>
#include <string>
#include <type_traits>

using FieldT = libff::Fr<libsnark::default_r1cs_ppzksnark_pp>;
using Hash = Mimc512F;
using GadHash = mimc512f_two_to_one_hash_gadget;
#else
template<typename FieldT, typename Hash, typename GadHash>
#endif
class PABR_Gadget : public libsnark::gadget<FieldT>
{
public:
    using super = libsnark::gadget<FieldT>;
    using PbVar = libsnark::pb_variable<FieldT>;
    using DigVar = std::pair<PbVar, PbVar>;
    using Protoboard = libsnark::protoboard<FieldT>;
    using Constraint = libsnark::r1cs_constraint<FieldT>;
    using LC = libsnark::linear_combination<FieldT>;

    static constexpr size_t DIGEST_BITS = Hash::DIGEST_SIZE * CHAR_BIT;

private:
    DigVar trans;
    DigVar other;
    std::vector<DigVar> middle;
    std::vector<DigVar> otherx;
    std::vector<DigVar> otherxm;
    std::vector<DigVar> inter;
    std::vector<DigVar> interx;
    std::vector<DigVar> interxm;
    std::vector<GadHash> foo_hash;
    size_t height;
    size_t leaves_n;
    size_t idx;

public:
    const DigVar out;

    PABR_Gadget(Protoboard &pb,                                                       //
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
        leaves_n{1ULL << (height - 1)}, //
        idx{trans_idx},                 //
        out{out}
    {
        // our node is a leaf
        if (trans_idx < leaves_n)
        {
            interx.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, other, trans, interx[0], FMT(""));
            else
                foo_hash.emplace_back(pb, trans, other, interx[0], FMT(""));

            interxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            otherxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});

            trans_idx >>= 1;
            inter.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, otherxm[0], interxm[0], inter[0], FMT(""));
            else
                foo_hash.emplace_back(pb, interxm[0], otherxm[0], inter[0], FMT(""));
        }
        // our node is a middle node
        else
        {
            // translate trans_idx to be used as a path along the tree
            trans_idx >>= 1;
            do
            {
                trans_idx -= leaves_n;
                leaves_n /= 2;
                --height;
            } while (trans_idx >= leaves_n);

            interxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            otherxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});

            inter.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            foo_hash.emplace_back(pb, interxm[0], otherxm[0], inter[0], FMT(""));
        }

        for (size_t i = 1; i < height - 2; ++i)
        {
            interx.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            interxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            otherxm.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});
            inter.emplace_back(PbVar{pb, FMT("")}, PbVar{pb, FMT("")});

            trans_idx >>= 1;
            if (trans_idx & 1)
                foo_hash.emplace_back(pb, otherxm[i], interx[i], inter[i], FMT(""));
            else
                foo_hash.emplace_back(pb, interx[i], otherxm[i], inter[i], FMT(""));
        }
    }

    inline size_t constrain(const LC &a, const LC &b, const LC &c)
    {
        this->pb.add_r1cs_constraint(Constraint(a, b, c), FMT(""));

        return 1;
    }

    void generate_r1cs_constraints()
    {
        size_t trans_idx = idx;

        // our node is a leaf
        if (trans_idx < leaves_n)
        {
            constrain(interx[0].first + middle[0].first, 1, interxm[0].first);
            constrain(interx[0].second + middle[0].second, 1, interxm[0].second);

            constrain(otherx[0].first + middle[0].second, 1, otherxm[0].first);
            constrain(otherx[0].second + middle[0].second, 1, otherxm[0].second);
            trans_idx >>= 1;
        }
        // our node is a middle node
        else
        {
            // translate trans_idx to be used as a path along the tree
            trans_idx >>= 1;
            do
            {
                trans_idx -= leaves_n;
                leaves_n /= 2;
                --height;
            } while (trans_idx >= leaves_n);

            constrain(trans.first + other.first, 1, interxm[0].first);
            constrain(trans.second + other.second, 1, interxm[0].second);

            constrain(trans.first + otherx[0].first, 1, otherxm[0].first);
            constrain(trans.second + otherx[0].second, 1, otherxm[0].second);
        }

        for (size_t i = 1; i < height - 2; ++i)
        {
            // the node was a right leaf, the remainder is interx[0]
            if (trans_idx & 1)
            {
                constrain(inter[i - 1].first + interx[i - 1].first, 1, interx[i].first);
                constrain(inter[i - 1].second + interx[i - 1].second, 1, interx[i].second);
            }
            // the node was a left leaf or a middle, the remainder is otherx[0]
            else
            {
                constrain(inter[i - 1].first + otherx[i - 1].first, 1, interx[i].first);
                constrain(inter[i - 1].second + otherx[i - 1].second, 1, interx[i].second);
            }

            trans_idx >>= 1;
            constrain(interx[i].first + middle[i].first, 1, interxm[i].first);
            constrain(interx[i].second + middle[i].second, 1, interxm[i].second);

            constrain(otherx[i].first + middle[i].first, 1, otherxm[i].first);
            constrain(otherx[i].second + middle[i].second, 1, otherxm[i].second);
        }
        if (trans_idx & 1)
        {
            constrain(inter.back().first + interx.back().first, 1, out.first);
            constrain(inter.back().second + interx.back().second, 1, out.second);
        }
        else
        {
            constrain(inter.back().first + otherx.back().first, 1, out.first);
            constrain(inter.back().second + otherx.back().second, 1, out.second);
        }

        for (auto &&x : foo_hash)
            x.generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        size_t trans_idx = idx;
        FieldT remain;

        // our node is a leaf
        if (trans_idx < leaves_n)
        {
            foo_hash[0].generate_r1cs_witness();

            trans_idx >>= 1;
            if (trans_idx & 1)
                remain = pb.val(interx[0]);
        }
    }
};
