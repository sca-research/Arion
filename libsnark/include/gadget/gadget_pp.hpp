#pragma once

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

template<typename FieldT>
class GadgetPP : public libsnark::gadget<FieldT>
{
public:
    using super = libsnark::gadget<FieldT>;

    using super::super;

protected:
    using LC = libsnark::linear_combination<FieldT>;

    inline size_t constrain(const LC &x, const LC &y, const LC &z)
    {
        this->pb.add_r1cs_constraint(libsnark::r1cs_constraint<FieldT>(x, y, z), FMT(""));

        return 1;
    }

    inline FieldT val(const libsnark::pb_variable<FieldT> &x) const { return this->pb.val(x); }
    inline FieldT &val(const libsnark::pb_variable<FieldT> &x) { return this->pb.val(x); }

    inline FieldT val(const libsnark::pb_linear_combination<FieldT> &x) const
    {
        x.evaluate(this->pb);
        return this->pb.lc_val(x);
    }

    inline FieldT &val(const libsnark::pb_linear_combination<FieldT> &x)
    {
        x.evaluate(this->pb);
        return this->pb.lc_val(x);
    }
};
