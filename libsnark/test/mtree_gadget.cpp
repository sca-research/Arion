#include "gadget/mtree_gadget.hpp"
#include "gadget/griffin/griffin_gadget.hpp"
#include "gadget/mimc256/mimc256_gadget.hpp"
#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "gadget/mimc512f2k/mimc512f2k_gadget.hpp"
#include "gadget/poseidon5/poseidon5_gadget.hpp"
#include "gadget/sha256/sha256_gadget_pp.hpp"
#include "gadget/sha512/sha512_gadget_pp.hpp"
#include "gadget/arion/arion_gadget.hpp"
#include "util/measure.hpp"
#include "tree/mtree.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>
#include <omp.h>

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

template<typename GadTree, bool full_tree = false>
bool test_mtree(size_t trans_idx = 0)
{
    static constexpr size_t HEIGHT = GadTree::HEIGHT;

    using DigVar = typename GadTree::DigVar;
    using Level = typename GadTree::Level;
    using GadHash = typename GadTree::GadHash;
    using Hash = typename GadHash::Hash;
    using Tree = std::conditional_t<full_tree, MTree<HEIGHT, Hash>, MTreePath<HEIGHT, Hash>>;
    using Node = typename Tree::Node;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    // Build tree
    std::vector<uint8_t> data(Tree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Tree tree{data.begin(), data.end()};

    std::cout << '\n' << tree;

    // There is no need to extract as FieldVariable does that for us

    // Test Gadget
    libsnark::protoboard<FieldT> pb;

    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<Level> other;
    PbVariablePP<FieldT> idx{pb, FMT("idx")};

    for (size_t i = 0; i < HEIGHT - 1; ++i)
        other.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("other_%llu", i)));

    GadTree gadget{pb, out, trans, other, idx, FMT("merkle_tree")};

    pb.set_input_sizes(DIGEST_VARS);
    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        for (size_t j = 0; j < other[i].size(); ++j)
            other[i][j].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    //    out.generate_r1cs_witness(tree.digest());
    trans.generate_r1cs_witness(tree.get_node(trans_idx)->get_digest());
    pb.val(idx) = trans_idx;

    const Node *aux = tree.get_node(trans_idx)->get_f();
    for (size_t i = 0; i < other.size(); ++i, aux = aux->get_f())
        for (size_t j = 0; j < other[i].size(); ++j)
            other[i][j].generate_r1cs_witness(aux->get_c(j)->get_digest());

    gadget.generate_r1cs_witness();

    std::string vanilla_dump{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
    std::string zkp_dump;

    if constexpr (GadTree::HASH_ISBOOLEAN)
    {
        uint8_t buff[Hash::DIGEST_SIZE]{};
        std::vector<bool> buff_bv(DIGEST_VARS);

        for (size_t i = 0; i < DIGEST_VARS; ++i)
            buff_bv[i] = pb.val(out.bits[i]).as_ulong();
        pack_bits(buff, buff_bv);

        zkp_dump = hexdump(buff);
    }
    else
    {
        for (auto &&x : out)
            zkp_dump += hexdump(pb.val(x).as_bigint());
    }

    std::cout << "\nVanilla output:\t" << vanilla_dump << '\n';
    std::cout << "ZKP output:\t" << zkp_dump << '\n';

    bool result = vanilla_dump == zkp_dump;

    auto keypair{libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system())};
    auto proof{
        libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(), pb.auxiliary_input())};

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

static bool run_tests()
{
    static constexpr size_t TREE_HEIGHT = 4;

    bool check = true;
    bool all_check = true;
    std::cout << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    ppT::init_public_params();

    std::cout << "SHA256... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha256Gadget<FieldT>>, true>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA256... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha256Gadget<FieldT>>>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC256... ";
    std::cout.flush();
    {
        //      check = test_pmtree<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC256... ";
    std::cout.flush();
    {
        //     check = test_pmtree_path<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA512... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha512Gadget<FieldT>>>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F... ";
    std::cout.flush();
    {
        //      check = test_pmtree_path<TREE_HEIGHT, GadMimc512F>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F2K... ";
    std::cout.flush();
    {
        //      check = test_pmtree_path<TREE_HEIGHT, GadMimc512F2K>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Griffin... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, GriffinGadget<Griffin<FieldT, 2, 1>>>>();
    }
    std::cout << check << '\n';
    all_check &= check;


    std::cout << "Path Poseidon5... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Poseidon5Gadget<Poseidon5<FieldT, 2, 1>>>>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Arion... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, ArionGadget<Arion<FieldT, 2, 1>>>>();
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
#endif

    return 0;
}
