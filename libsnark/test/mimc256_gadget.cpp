#include "gadget/mimc256/mimc256_gadget.hpp"
#include "hash/mimc256.hpp"
#include "util/measure.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>


using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

bool test()
{
    using DigVar = FieldVariable<FieldT>;
    using GadHash = mimc256_two_to_one_hash_gadget<FieldT>;
    using Hash = GadHash::Hash;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    std::vector<uint8_t> block(Hash::BLOCK_SIZE);
    std::vector<uint8_t> digest(Hash::DIGEST_SIZE);
    std::generate(block.begin(), block.end(), std::ref(rng));
    std::string vanilla_dump;
    std::string zksnark_dump;

    Hash::hash_oneblock(digest.data(), block.data());

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar left{pb, DIGEST_VARS, FMT("trans")};
    DigVar right{pb, DIGEST_VARS, FMT("trans")};

    pb.set_input_sizes(DIGEST_VARS);

    GadHash gadget{pb, left, right, out, ""};

    out.generate_r1cs_constraints();
    left.generate_r1cs_constraints();
    right.generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();


    left.generate_r1cs_witness(block);
    right.generate_r1cs_witness(block.data() + Hash::DIGEST_SIZE, Hash::DIGEST_SIZE);
    gadget.generate_r1cs_witness();

    vanilla_dump = hexdump(digest.begin(), digest.end());

    for (auto &&x : out)
        zksnark_dump += hexdump(pb.val(x).as_bigint());

    std::cout << '\n' << "Vanilla output:\t" << vanilla_dump << '\n';
    std::cout << "ZKP output:\t" << zksnark_dump << '\n';

    bool result = vanilla_dump == zksnark_dump;
    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

static bool run_tests()
{
    bool check = true;
    bool all_check = true;
    std::cout << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    ppT::init_public_params();


    std::cout << "MiMC256... ";
    std::cout.flush();
    check = test();
    std::cout << check << '\n';
    all_check &= check;


    return all_check;
}

int main()
{
    std::cout << "\n==== Testing MiMC256 Gadget ====\n";

    bool all_check = false;
    measure([&]() { all_check = run_tests(); });

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
