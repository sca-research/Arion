#pragma once

#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>

#include "gadget/gadget_pp.hpp"

template<typename FieldT>
class XORGadget : public GadgetPP<FieldT>
{
public:
    using super = GadgetPP<FieldT>;
    using PbVar = libsnark::pb_variable<FieldT>;

    using super::constrain;
    using super::val;

    const PbVar x;
    const PbVar y;
    const PbVar out;

    XORGadget(libsnark::protoboard<FieldT> &pb, const PbVar &x, const PbVar &y, const PbVar &out,
              const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}
    {}

    void generate_r1cs_constraints()
    {
        // x(1 - x) = 0, i.e. x must be 0 or 1
        constrain(x, 1 - x, 0);
        constrain(y, 1 - y, 0);

        // (x + y) - out = (x + x)y  i.e. out = x ^ y
        constrain(x + x, y, x + y - out);
    }

    void generate_r1cs_witness() { val(out) = val(x) + val(y) - (val(x) + val(x)) * val(y); }
};

template<typename FieldT>
class LongXORGadget : public libsnark::gadget<FieldT>
{
public:
    using super = libsnark::gadget<FieldT>;
    using DigVar = libsnark::digest_variable<FieldT>;

    const DigVar x;
    const DigVar y;
    const DigVar out;

private:
    std::vector<XORGadget<FieldT>> xor_gad;

public:
    LongXORGadget(libsnark::protoboard<FieldT> &pb, const DigVar &x, const DigVar &y,
                  const DigVar &out, const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}
    {
        for (size_t i = 0; i < x.digest_size; ++i)
            xor_gad.emplace_back(pb, x.bits[i], y.bits[i], out.bits[i],
                                 FMT(annotation_prefix, "_xor"));
    }

    void generate_r1cs_constraints()
    {
        for (auto &&x : xor_gad)
            x.generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        for (auto &&x : xor_gad)
            x.generate_r1cs_witness();
    }
};
