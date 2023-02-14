#pragma once

#include "gadget/digest_variable_pp.hpp"
#include "gadget/sha512/sha512_gadget.hpp"
#include "hash/sha512.hpp"
#include <array>

template<typename FieldT>
class Sha512Gadget : public libsnark::sha512::sha512_two_to_one_hash_gadget<FieldT>
{
public:
    using super = libsnark::sha512::sha512_two_to_one_hash_gadget<FieldT>;
    using Hash = Sha512;
    using Field = FieldT;
    using DigVar = DigestVariablePP<Field>;
    using BlockVar = std::array<DigVar, 2>;

    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = DIGEST_SIZE * CHAR_BIT;

    Sha512Gadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                 const std::string &ap) :
        super(pb, in[0], in[1], out, ap)
    {}

    void generate_r1cs_constraints() { super::generate_r1cs_constraints(true); }
    void generate_r1cs_witness() { super::generate_r1cs_witness(); }
};

template<typename FieldT>
class sha512_two_to_one_hash_gadget : public Sha512Gadget<FieldT>
{
public:
    using super = Sha512Gadget<FieldT>;
    using Hash = typename super::Hash;
    using DigVar = typename super::DigVar;
    using Field = typename super::Field;
    using BlockVar = typename super::BlockVar;

    sha512_two_to_one_hash_gadget(libsnark::protoboard<Field> &pb, const DigVar &x, const DigVar &y,
                                  const DigVar &out, const std::string &ap) :
        super{pb, BlockVar{x, y}, out, ap}
    {}
};
