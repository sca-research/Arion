#pragma once

#include <libsnark/gadgetlib1/pb_variable.hpp>

template<typename FieldT>
class PbVariablePP : public libsnark::pb_variable<FieldT>
{
public:
    using super = libsnark::pb_variable<FieldT>;
    using super::super;

    PbVariablePP(libsnark::protoboard<FieldT> &pb, const std::string &ap) : super{}
    {
        this->allocate(pb, FMT(ap));
    };
};

template<typename FieldT>
class PbLCPP : public libsnark::pb_linear_combination<FieldT>
{
public:
    using super = libsnark::pb_linear_combination<FieldT>;
    using super::super;

    PbLCPP(libsnark::protoboard<FieldT> &pb, const libsnark::linear_combination<FieldT> &lc) :
        super{}
    {
        this->assign(pb, lc);
    };
};
