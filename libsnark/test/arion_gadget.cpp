#include "gadget/arion/arion_gadget.hpp"
#include "util/array_utils.hpp"
#include "hash/arion.hpp"
#include "util/measure.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>


using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

bool test()
{
    using GadHash = ArionGadget<Arion<FieldT, 7, 3>>;
    using Hash = GadHash::Hash;
    using DigVar = GadHash::DigVar;
    using BlockVar = GadHash::BlockVar;

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
    BlockVar in{make_uniform_array<BlockVar>(pb, DIGEST_VARS, FMT("trans"))};

    pb.set_input_sizes(DIGEST_VARS);

    GadHash gadget{pb, in, out, ""};

    out.generate_r1cs_constraints();

    for (size_t i = 0; i < in.size(); ++i)
        in[i].generate_r1cs_constraints();

    gadget.generate_r1cs_constraints();

    for (size_t i = 0; i < in.size(); ++i)
        in[i].generate_r1cs_witness(block.data() + i * Hash::DIGEST_SIZE, Hash::DIGEST_SIZE);
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


    std::cout << "Hashing... ";
    std::cout.flush();
    check = test();
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Arion Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
