#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"

template<typename FieldT>
class AddGadget : public GadgetPP<FieldT>
{
public:
    using super = GadgetPP<FieldT>;
    using PbVar = libsnark::pb_variable<FieldT>;

    using super::constrain;
    using super::val;

    const PbVar x;
    const PbVar y;
    const PbVar out;

    AddGadget(libsnark::protoboard<FieldT> &pb, const PbVar &x, const PbVar &y, const PbVar &out,
              const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}
    {}

    void generate_r1cs_constraints() { constrain(x + y, 1, out); }

    void generate_r1cs_witness() { val(out) = val(x) + val(y); }
};

template<typename FieldT>
class LongAddGadget : public GadgetPP<FieldT>
{
public:
    using super = GadgetPP<FieldT>;
    using DigVar = FieldVariable<FieldT>;

    const DigVar x;
    const DigVar y;
    const DigVar out;

private:
    std::vector<AddGadget<FieldT>> add_gad;

public:
    LongAddGadget(libsnark::protoboard<FieldT> &pb, const DigVar &x, const DigVar &y,
                  const DigVar &out, const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}
    {
        for (size_t i = 0; i < x.size(); ++i)
            add_gad.emplace_back(pb, x[i], y[i], out[i], FMT(annotation_prefix, "_add"));
    }

    void generate_r1cs_constraints()
    {
        for (auto &&x : add_gad)
            x.generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        for (auto &&x : add_gad)
            x.generate_r1cs_witness();
    }
};
