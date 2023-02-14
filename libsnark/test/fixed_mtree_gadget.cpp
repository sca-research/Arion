#include "gadget/fixed_mtree_gadget.hpp"
#include "gadget/griffin/griffin_gadget.hpp"
#include "gadget/mimc256/mimc256_gadget.hpp"
#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "gadget/mimc512f2k/mimc512f2k_gadget.hpp"
#include "gadget/poseidon5/poseidon5_gadget.hpp"
#include "gadget/sha256/sha256_gadget_pp.hpp"
#include "gadget/sha512/sha512_gadget.hpp"
#include "gadget/arion/arion_gadget.hpp"
#include "tree/fixed_mtree.hpp"
#include "util/measure.hpp"

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
using GadSha512 = libsnark::sha512::sha512_two_to_one_hash_gadget<FieldT>;

using GadGriffin = griffin_two_to_one_hash_gadget<FieldT>;
using GadMimc256 = mimc256_two_to_one_hash_gadget<FieldT>;
using GadMimc512F = mimc512f_two_to_one_hash_gadget<FieldT>;
using GadMimc512F2K = mimc512f2k_two_to_one_hash_gadget<FieldT>;
using GadPoseidon5 = poseidon5_two_to_one_hash_gadget<FieldT>;
using GadArion = arion_two_to_one_hash_gadget<FieldT>;

template<size_t tree_height, typename GadHash>
bool test_mtree()
{
    using Hash = typename GadHash::Hash;
    using Mtree = FixedMTree<tree_height, Hash>;
    using DigVar = libsnark::digest_variable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    // Build tree
    std::vector<uint8_t> data(Mtree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Mtree tree{data.begin(), data.end()};

    std::cout << '\n' << tree << '\n';

    // Extract our transaction
    libff::bit_vector trans_bv(DIGEST_VARS);
    unpack_bits(trans_bv, tree.get_node(TRANS_IDX)->get_digest());

    // Extract other transaction/middle nodes
    std::vector<libff::bit_vector> other_bv(tree_height - 1, libff::bit_vector(DIGEST_VARS));
    for (size_t i = 0, j = 1; i < other_bv.size(); ++i, j += 1ULL << (tree_height - i))
        unpack_bits(other_bv[i], tree.get_node(j)->get_digest());

    // Extract output node
    libff::bit_vector out_bv(DIGEST_VARS);
    unpack_bits(out_bv, tree.digest());

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<DigVar> other;

    for (size_t i = 0; i < tree_height - 1; ++i)
        other.emplace_back(pb, DIGEST_VARS, FMT("other_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    FixedMTreeGadget<GadHash> gadget{pb, out, trans, other, TRANS_IDX, FMT("merkle_tree")};

    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();


    //    out.generate_r1cs_witness(out_bv);
    trans.generate_r1cs_witness(trans_bv);
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_witness(other_bv[i]);
    gadget.generate_r1cs_witness();

    bool result;
    {
        uint8_t buff[Hash::DIGEST_SIZE]{};
        std::vector<bool> buff_bv(DIGEST_VARS);

        for (size_t i = 0; i < DIGEST_VARS; ++i)
            buff_bv[i] = pb.val(out.bits[i]).as_ulong();
        pack_bits(buff, buff_bv);

        std::string vanilla_dump{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
        std::string zkp_dump{hexdump(buff)};

        //std::cout << "\nVanilla output:\t" << vanilla_dump << '\n';
        //std::cout << "ZKP output:\t" << zkp_dump << '\n';

        result = vanilla_dump == zkp_dump;
    }

    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

template<size_t tree_height, typename GadHash>
bool test_mtree_path()
{
    using Hash = typename GadHash::Hash;
    using Mtree = FixedMTreePath<tree_height, Hash>;
    using DigVar = libsnark::digest_variable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    // Build tree
    std::vector<uint8_t> data(Mtree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Mtree tree{data.begin(), data.end()};

    std::cout << '\n' << tree << '\n';

    // Extract our transaction
    libff::bit_vector trans_bv(DIGEST_VARS);
    unpack_bits(trans_bv, data.data());

    // Extract other transaction/middle nodes
    std::vector<libff::bit_vector> other_bv(tree_height - 1, libff::bit_vector(DIGEST_VARS));
    for (size_t i = 0; i < other_bv.size(); ++i)
        unpack_bits(other_bv[i], data.data() + (i + 1) * Hash::DIGEST_SIZE);

    // Extract output node
    libff::bit_vector out_bv(DIGEST_VARS);
    unpack_bits(out_bv, tree.digest());

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<DigVar> other;

    for (size_t i = 0; i < tree_height - 1; ++i)
        other.emplace_back(pb, DIGEST_VARS, FMT("other_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    FixedMTreeGadget<GadHash> gadget{pb, out, trans, other, TRANS_IDX, FMT("merkle_tree")};

    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();


    //    out.generate_r1cs_witness(out_bv);
    trans.generate_r1cs_witness(trans_bv);
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_witness(other_bv[i]);
    gadget.generate_r1cs_witness();

    bool result;
    {
        uint8_t buff[Hash::DIGEST_SIZE]{};
        std::vector<bool> buff_bv(DIGEST_VARS);

        for (size_t i = 0; i < DIGEST_VARS; ++i)
            buff_bv[i] = pb.val(out.bits[i]).as_ulong();
        pack_bits(buff, buff_bv);

        std::string vanilla_dump{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
        std::string zkp_dump{hexdump(buff)};

        //        std::cout << "\nVanilla output:\t" << vanilla_dump << '\n';
        //        std::cout << "ZKP output:\t" << zkp_dump << '\n';

        result = vanilla_dump == zkp_dump;
    }

    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

template<size_t tree_height, typename GadHash>
bool test_pmtree()
{
    using Hash = typename GadHash::Hash;
    using Mtree = FixedMTree<tree_height, Hash>;
    using DigVar = FieldVariable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{4};

    // Build tree
    std::vector<uint8_t> data(Mtree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Mtree tree{data.begin(), data.end()};

    std::cout << '\n' << tree;

    // There is no need to extract as FieldVariable does that for us

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<DigVar> other;

    for (size_t i = 0; i < tree_height - 1; ++i)
        other.emplace_back(pb, DIGEST_VARS, FMT("other_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    FixedMTreeGadget<GadHash> gadget{pb, out, trans, other, TRANS_IDX, FMT("merkle_tree")};

    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    //    out.generate_r1cs_witness(tree.digest());
    trans.generate_r1cs_witness(tree.get_node(TRANS_IDX)->get_digest());
    for (size_t i = 0, j = 1; i < other.size(); ++i, j += 1ULL << (tree_height - i))
        other[i].generate_r1cs_witness(tree.get_node(j)->get_digest());
    gadget.generate_r1cs_witness();


    bool result;
    {
        std::string vanilla{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
        std::string zkp;

        for (auto &&x : out)
            zkp += hexdump(pb.val(x).as_bigint());

        //        std::cout << "\nVanilla output:\t" << vanilla << '\n';
        //        std::cout << "ZKP output:\t" << zkp << '\n';

        result = vanilla == zkp;
    }

    auto keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system());
    auto proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                      pb.auxiliary_input());

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

template<size_t tree_height, typename GadHash>
bool test_pmtree_path()
{
    using Hash = typename GadHash::Hash;
    using Mtree = FixedMTreePath<tree_height, Hash>;
    using DigVar = FieldVariable<FieldT>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    // Build tree
    std::vector<uint8_t> data(Mtree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Mtree tree{data.begin(), data.end()};

    std::cout << '\n' << tree << '\n';

    // There is no need to extract anything as FieldVariable does that for us

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<DigVar> other;

    for (size_t i = 0; i < tree_height - 1; ++i)
        other.emplace_back(pb, DIGEST_VARS, FMT("other_%llu", i));

    pb.set_input_sizes(DIGEST_VARS);
    FixedMTreeGadget<GadHash> gadget{pb, out, trans, other, 0, FMT("merkle_tree")};

    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    out.generate_r1cs_witness(tree.digest());
    trans.generate_r1cs_witness(data.data());
    for (size_t i = 0; i < other.size(); ++i)
        other[i].generate_r1cs_witness(data.data() + (i + 1) * Hash::DIGEST_SIZE);
    gadget.generate_r1cs_witness();

    bool result;
    {
        std::string vanilla{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
        std::string zkp;

        for (auto &&x : out)
            zkp += hexdump(pb.val(x).as_bigint());

        //       std::cout << "Vanilla output:\t" << hexdump(tree.digest(), Hash::DIGEST_SIZE) << '\n';
        //       std::cout << "ZKP output:\t" << zkp << '\n';

        result = vanilla == zkp;
    }

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

    std::cout << "SHA256... ";
    std::cout.flush();
    {
        check = test_mtree<TREE_HEIGHT, GadSha256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA256... ";
    std::cout.flush();
    {
        check = test_mtree_path<TREE_HEIGHT, GadSha256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC256... ";
    std::cout.flush();
    {
        check = test_pmtree<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC256... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA512... ";
    std::cout.flush();
    {
        check = test_mtree_path<TREE_HEIGHT, GadSha512>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadMimc512F>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F2K... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadMimc512F2K>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Griffin... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadGriffin>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Poseidon5... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadPoseidon5>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Arion... ";
    std::cout.flush();
    {
        check = test_pmtree_path<TREE_HEIGHT, GadArion>();
    }
    std::cout << check << '\n';
    all_check &= check;


    return all_check;
}

int main()
{
    std::cout << "\n==== Testing MerkleTree Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
    std::cout << "Measuring SHA256 performance:\n";
    measure([]() { test_mtree<TREE_HEIGHT, GadSha256>(); });
    std::cout << "Measuring SHA512 performance:\n";
    measure([]() { test_mtree<TREE_HEIGHT, GadSha512>(); });
#endif

    return 0;
}
