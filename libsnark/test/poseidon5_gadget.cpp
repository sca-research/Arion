#include "gadget/poseidon5/poseidon5_gadget.hpp"
#include "util/array_utils.hpp"
#include "hash/poseidon5.hpp"
#include "util/measure.hpp"
#include <fstream>
#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

bool test()
{
    using GadHash = Poseidon5Gadget<Poseidon5<FieldT, 2, 1>>;
    using DigVar = GadHash::DigVar;
    using BlockVar = GadHash::BlockVar;
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
    BlockVar in{make_uniform_array<BlockVar>(pb, DIGEST_VARS, FMT("trans"))};
    DigVar out{pb, DIGEST_VARS, FMT("out")};

    pb.set_input_sizes(DIGEST_VARS);

    GadHash gadget{pb, in, out, ""};

    out.generate_r1cs_constraints();
    for (auto &&x : in)
        x.generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    for (size_t i = 0; i < Hash::RATE; ++i)
        in[i].generate_r1cs_witness(block.data() + i * Hash::DIGEST_SIZE, Hash::DIGEST_SIZE);
    gadget.generate_r1cs_witness();

    vanilla_dump = hexdump(digest);
    for (auto &&x : out)
        zksnark_dump += hexdump(pb.val(x));
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


    std::cout << "Hashing... ";
    std::cout.flush();
    check = test();
    std::cout << check << '\n';
    all_check &= check;


    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Poseidon5 Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
