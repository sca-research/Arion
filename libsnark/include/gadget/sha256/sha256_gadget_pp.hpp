#pragma once

#include "hash/sha256.hpp"
#include <array>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>
#include "gadget/digest_variable_pp.hpp"

template<typename FieldT>
class Sha256Gadget : public libsnark::sha256_two_to_one_hash_gadget<FieldT>
{
public:
    using super = libsnark::sha256_two_to_one_hash_gadget<FieldT>;
    using Hash = Sha256;
    using Field = FieldT;
    using DigVar = DigestVariablePP<Field>;
    using BlockVar = std::array<DigVar, 2>;

    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = DIGEST_SIZE * CHAR_BIT;

    Sha256Gadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                 const std::string &ap) :
        super(pb, in[0], in[1], out, ap)
    {}

    void generate_r1cs_constraints() { super::generate_r1cs_constraints(true); }
    void generate_r1cs_witness() { super::generate_r1cs_witness(); }
};

template<typename FieldT>
class sha256_two_to_one_hash_gadget : public Sha256Gadget<FieldT>
{
public:
    using super = Sha256Gadget<FieldT>;
    using Hash = typename super::Hash;
    using Field = typename super::Field;
    using DigVar = typename super::DigVar;
    using BlockVar = typename super::BlockVar;

    sha256_two_to_one_hash_gadget(libsnark::protoboard<Field> &pb, const DigVar &x, const DigVar &y,
                                  const DigVar &out, const std::string &ap) :
        super{pb, BlockVar{x, y}, out, ap}
    {}

    using super::generate_r1cs_constraints;
    using super::generate_r1cs_witness;
};
