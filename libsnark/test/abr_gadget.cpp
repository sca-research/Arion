#include "gadget/abr_gadget.hpp"
#include "gadget/mimc256/mimc256_gadget.hpp"
#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "gadget/mimc512f2k/mimc512f2k_gadget.hpp"
#include "gadget/sha256/sha256_gadget_pp.hpp"
#include "gadget/sha512/sha512_gadget_pp.hpp"
#include "tree/fixed_abr.hpp"
#include "util/measure.hpp"
#include "hash/sha256.hpp"
#include "hash/sha512.hpp"


#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>
#include <omp.h>

static constexpr size_t TRANS_IDX = 0;
static constexpr size_t TREE_HEIGHT = 4;

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

using GadSha256 = sha256_two_to_one_hash_gadget<FieldT>;
using GadSha512 = sha512_two_to_one_hash_gadget<FieldT>;
using GadMimc256 = mimc256_two_to_one_hash_gadget<FieldT>;
using GadMimc512F = mimc512f_two_to_one_hash_gadget<FieldT>;
using GadMimc512F2K = mimc512f2k_two_to_one_hash_gadget<FieldT>;

template<size_t tree_height, typename GadHash>
bool test_tRee()
{
    using Hash = typename GadHash::Hash;
    using Abr = FixedAbr<tree_height, Hash>;
    using DigVar = libsnark::digest_variable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};


    // Build tree
    std::vector<uint8_t> data(Abr::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Abr tree{data.begin(), data.end()};

    std::cout << '\n' << tree << '\n';

    // Extract our transaction
    libff::bit_vector trans_bv(DIGEST_VARS);
    unpack_bits(trans_bv, tree.get_node(TRANS_IDX)->get_digest());
    std::cout << "trans_bv: " << hexdump(trans_bv) << '\n';

    // Extract other transaction node
    libff::bit_vector other_bv(DIGEST_VARS);
    unpack_bits(other_bv, tree.get_node(TRANS_IDX + 1)->get_digest());

    std::cout << "other_bv: " << hexdump(other_bv) << '\n';

    // Extract middle nodes
    std::vector<libff::bit_vector> middle_bv(tree_height - 2, libff::bit_vector(DIGEST_VARS));
    for (size_t i = 0, j = Abr::LEAVES_N; i < middle_bv.size();
         ++i, j += 1ULL << (tree_height - 2 - i))
        unpack_bits(middle_bv[i], tree.get_node(j)->get_digest());

    for (auto &&x : middle_bv)
        std::cout << "middle_bv: " << hexdump(x) << '\n';

    // Extract otherx nodes
    std::vector<libff::bit_vector> otherx_bv(tree_height - 2, libff::bit_vector(DIGEST_VARS));
    for (size_t i = 0, j = Abr::INPUT_N + 1; i < otherx_bv.size();
         ++i, j += 1ULL << (tree_height - 1 - i))
        unpack_bits(otherx_bv[i], tree.get_node(j)->get_digest());

    for (auto &&x : otherx_bv)
        std::cout << "otherx_bv: " << hexdump(x) << '\n';

    // Extract output node
    libff::bit_vector out_bv(DIGEST_VARS);
    unpack_bits(out_bv, tree.digest());

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    DigVar other{pb, DIGEST_VARS, FMT("trans")};

    std::vector<DigVar> middle;
    for (size_t i = 0; i < middle_bv.size(); ++i)
        middle.emplace_back(pb, DIGEST_VARS, FMT("middle_%llu", i));

    std::vector<DigVar> otherx;
    for (size_t i = 0; i < otherx_bv.size(); ++i)
        otherx.emplace_back(pb, DIGEST_VARS, FMT("otherx_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    ABR_Gadget<FieldT, GadHash> gadget{
        pb, out, trans, other, middle, otherx, TRANS_IDX, tree_height, FMT("merkle_tree")};

    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    other.generate_r1cs_constraints();
    for (auto &&x : middle)
        x.generate_r1cs_constraints();
    for (auto &&x : otherx)
        x.generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    out.generate_r1cs_witness(out_bv);

    trans.generate_r1cs_witness(trans_bv);
    other.generate_r1cs_witness(other_bv);

    for (size_t i = 0; i < middle.size(); ++i)
        middle[i].generate_r1cs_witness(middle_bv[i]);

    for (size_t i = 0; i < otherx.size(); ++i)
        otherx[i].generate_r1cs_witness(otherx_bv[i]);

    gadget.generate_r1cs_witness();

    {
        std::cout << '\n'
                  << "Vanilla output:\t" << hexdump(tree.digest(), Hash::DIGEST_SIZE) << '\n';

        uint8_t buff[Hash::DIGEST_SIZE]{};
        std::vector<bool> buff_bv(DIGEST_VARS);

        for (size_t i = 0; i < DIGEST_VARS; ++i)
            buff_bv[i] = pb.val(out.bits[i]).as_ulong();
        pack_bits(buff, buff_bv);

        std::cout << "ZKP output:\t" << hexdump(buff) << '\n';
    }

    bool result;
    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                              proof);

    return result;
}

template<size_t tree_height, typename GadHash>
bool test_ptRee()
{
    using Hash = typename GadHash::Hash;
    using Abr = FixedAbr<tree_height, Hash>;
    using DigVar = FieldVariable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};


    // Build tree
    std::vector<uint8_t> data(Abr::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Abr tree{data.begin(), data.end()};

    // There is no need to extract anything as FieldVariable does that for us

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    DigVar other{pb, DIGEST_VARS, FMT("trans")};

    std::vector<DigVar> middle;
    for (size_t i = 0; i < tree_height - 2; ++i)
        middle.emplace_back(pb, DIGEST_VARS, FMT("middle_%llu", i));

    std::vector<DigVar> otherx;
    for (size_t i = 0; i < tree_height - 2; ++i)
        otherx.emplace_back(pb, DIGEST_VARS, FMT("otherx_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    ABR_Gadget<FieldT, GadHash> gadget{
        pb, out, trans, other, middle, otherx, TRANS_IDX, tree_height, FMT("merkle_tree")};

    std::cout << '\n' << tree << '\n';
    gadget.generate_r1cs_constraints();

    out.generate_r1cs_witness(tree.digest(), Hash::DIGEST_SIZE);

    trans.generate_r1cs_witness(tree.get_node(TRANS_IDX)->get_digest(), Hash::DIGEST_SIZE);
    std::cout << "trans: " << hexdump(pb.val(trans[0]).as_bigint()) << '\n';

    other.generate_r1cs_witness(tree.get_node(TRANS_IDX + 1)->get_digest(), Hash::DIGEST_SIZE);
    std::cout << "other: " << hexdump(pb.val(other[0]).as_bigint()) << '\n';


    for (size_t i = 0, j = Abr::LEAVES_N; i < middle.size();
         ++i, j += 1ULL << (tree_height - 2 - i))
        middle[i].generate_r1cs_witness(tree.get_node(j)->get_digest(), Hash::DIGEST_SIZE);

    for (auto &&x : middle)
    {
        std::cout << "middle: " << hexdump(pb.val(x[0]).as_bigint()) << '\n';
    }

    for (size_t i = 0, j = Abr::INPUT_N + 1; i < otherx.size();
         ++i, j += 1ULL << (tree_height - 1 - i))
        otherx[i].generate_r1cs_witness(tree.get_node(j)->get_digest(), Hash::DIGEST_SIZE);

    for (auto &&x : otherx)
        std::cout << "otherx: " << hexdump(pb.val(x[0]).as_bigint()) << '\n';

    gadget.generate_r1cs_witness();

    {
        std::cout << '\n'
                  << "Vanilla output:\t" << hexdump(tree.digest(), Hash::DIGEST_SIZE) << '\n';

        std::cout << "ZKP output:\t" << hexdump(pb.val(out[0]).as_bigint()) << '\n';
    }

    bool result;
    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
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

    std::cout << "SHA256... ";
    std::cout.flush();
    {
        check = test_tRee<TREE_HEIGHT,GadSha256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "SHA512... ";
    std::cout.flush();
    {
        check = test_tRee<TREE_HEIGHT, GadSha512>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC256... ";
    std::cout.flush();
    {
        check = test_ptRee<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC512F... ";
    std::cout.flush();
    {
        check = test_ptRee<TREE_HEIGHT, GadMimc512F>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC512f2k... ";
    std::cout.flush();
    {
        check = test_ptRee<TREE_HEIGHT, GadMimc512F2K>();
    }
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing ABR Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
