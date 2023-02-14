#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "hash/mimc512f.hpp"
#include "util/measure.hpp"

#include <fstream>
#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

bool test()
{
    using DigVar = FieldVariable<FieldT>;
    using GadHash = mimc512f_two_to_one_hash_gadget<FieldT>;
    using Hash = GadHash::Hash;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    std::array<uint8_t, Hash::BLOCK_SIZE> block{};
    //    std::array<uint8_t, Hash::BLOCK_SIZE> block =
    //        0x2c98b9b6ddab8feb3e5bc1cc829c4c866515b75863c0e867d422b452defbb85524cb8b07eb12d306b845d05d961c65516f5948c008c61f30003199790cee3be5039ca7efb1d40d8fadc0181f32445378bba01a2f4ab3fe16026d59b10d972ed70bf837b11b2867b5daab545fb90442ec210cb4e8d9c7d273fdb744017d25cc87_x;

    std::vector<uint8_t> digest(Hash::DIGEST_SIZE);
    std::generate(block.begin(), block.end(), std::ref(rng));
    std::string vanilla_dump;
    std::string zksnark_dump;

    // Initialize Gadgets
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar left{pb, DIGEST_VARS, FMT("left")};
    DigVar right{pb, DIGEST_VARS, FMT("right")};

    pb.set_input_sizes(DIGEST_VARS);

    GadHash gadget{pb, left, right, out, FMT("")};

    left.generate_r1cs_constraints();
    right.generate_r1cs_constraints();
    out.generate_r1cs_constraints();

    std::cout << '\n' << "Input:\t" << hexdump(block.begin(), block.end()) << '\n';
    Hash::hash_oneblock(digest.data(), block.data());
    vanilla_dump = hexdump(digest.begin(), digest.end());
    std::cout << "Vanilla:\t" << vanilla_dump << '\n';

    // Compute first level (SNARK)
    left.generate_r1cs_witness(block);
    right.generate_r1cs_witness(block.data() + Hash::DIGEST_SIZE, Hash::DIGEST_SIZE);
    gadget.generate_r1cs_witness();
    for (auto &&x : out)
        zksnark_dump += hexdump(pb.val(x).as_bigint());
    std::cout << "ZKSNARK:\t" << zksnark_dump << '\n';

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

    std::cout << "MiMC512F... ";
    std::cout.flush();
    check = test();
    std::cout << check << '\n';
    all_check &= check;


    return all_check;
}

int main()
{
    std::cout << "\n==== Testing MiMC512-Feistel Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
