#pragma once

#include "util/bit_pack.hpp"
#include <libsnark/gadgetlib1/gadgets/hashes/hash_io.hpp>


template<typename FieldT>
class DigestVariablePP : public libsnark::digest_variable<FieldT>
{
public:
    using super = libsnark::digest_variable<FieldT>;

    using super::super;

    template<typename... Args>
    void generate_r1cs_witness(Args &&...args)
    {
        super::generate_r1cs_witness(unpack_bits(std::forward<Args>(args)...));
    }

    const auto &operator[](size_t i) const { return this->bits[i]; }

    auto &operator[](size_t i) { return this->bits[i]; }

    auto begin() { return this->bits.begin(); }
    auto end() { return this->bits.end(); }
    auto begin() const { return this->bits.begin(); }
    auto end() const { return this->bits.end(); }
};
