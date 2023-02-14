#pragma once

#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

template<typename ppT>
class r1cs_ppzksnark_keypair
{
public:
    libsnark::r1cs_ppzksnark_proving_key<ppT> pk;
    libsnark::r1cs_ppzksnark_verification_key<ppT> vk;

    r1cs_ppzksnark_keypair() = default;
    r1cs_ppzksnark_keypair(const r1cs_ppzksnark_keypair<ppT> &other) = default;
    r1cs_ppzksnark_keypair(libsnark::r1cs_ppzksnark_proving_key<ppT> &&pk,
                           libsnark::r1cs_ppzksnark_verification_key<ppT> &&vk) :
        pk(std::move(pk)),
        vk(std::move(vk))
    {}
    r1cs_ppzksnark_keypair(const libsnark::r1cs_ppzksnark_keypair<ppT> &other) :
        pk(other.pk), vk(other.vk)
    {}

    r1cs_ppzksnark_keypair(libsnark::r1cs_ppzksnark_keypair<ppT> &&other) :
        pk(std::move(other.pk)), vk(std::move(other.vk))
    {}


    r1cs_ppzksnark_keypair(r1cs_ppzksnark_keypair<ppT> &&other) = default;

    r1cs_ppzksnark_keypair<ppT> &operator=(const r1cs_ppzksnark_keypair<ppT> &other) = default;
    r1cs_ppzksnark_keypair<ppT> &operator=(r1cs_ppzksnark_keypair<ppT> &&other) = default;
};
