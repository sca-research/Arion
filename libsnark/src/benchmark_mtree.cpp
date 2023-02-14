#include "gadget/griffin/griffin_gadget.hpp"
#include "gadget/mimc256/mimc256_gadget.hpp"
#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "gadget/mimc512f2k/mimc512f2k_gadget.hpp"
#include "gadget/mtree_gadget.hpp"
#include "gadget/poseidon5/poseidon5_gadget.hpp"
#include "gadget/sha256/sha256_gadget_pp.hpp"
#include "gadget/sha512/sha512_gadget_pp.hpp"
#include "gadget/arion/arion_gadget.hpp"
#include "gadget/arion_v2/arion_v2_gadget.hpp"
#include "r1cs/r1cs_ppzksnark_pp.hpp"
#include "tree/mtree.hpp"
#include "util/measure.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>
#include <omp.h>
#include <string>


static constexpr size_t MIN_HEIGHT = 4;
static constexpr size_t MAX_HEIGHT = 32 + 1; // the +1 is to highlight that the bound is exclusive

namespace fs = std::filesystem;

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

std::ofstream log_file;

template<size_t height, typename GadHash>
bool test_mtree(size_t trans_idx = 0)
{
    static constexpr size_t HEIGHT = height;

    using GadTree = MTreeGadget<height, GadHash>;
    using DigVar = typename GadTree::DigVar;
    using Level = typename GadTree::Level;
    using Hash = typename GadHash::Hash;
    using Tree = MTreePath<HEIGHT, Hash>;
    using Node = typename Tree::Node;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    double elap = 0;
    std::unique_ptr<Tree> tree;

    // Build tree
    elap = measure(
        [&]()
        {
            std::vector<uint8_t> data(Tree::INPUT_SIZE);
            std::generate(data.begin(), data.end(), std::ref(rng));
            tree = std::make_unique<Tree>(data.begin(), data.end());
        },
        1, 1, "Tree Generation", false);
    log_file << elap << '\t';
    log_file.flush();

    // Test Gadget
    libsnark::protoboard<FieldT> pb;

    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<Level> other;
    PbVariablePP<FieldT> idx{pb, FMT("idx")};
    std::unique_ptr<GadTree> gadget;

    for (size_t i = 0; i < HEIGHT - 1; ++i)
        other.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("other_%llu", i)));

    // Gadget construction
    elap = measure(
        [&]()
        { gadget = std::make_unique<GadTree>(pb, out, trans, other, idx, FMT("merkle_tree")); },
        1, 1, "Gadget construction", false);
    log_file << elap << '\t';
    log_file.flush();

    pb.set_input_sizes(DIGEST_VARS);

    // Constraint generation
    elap = measure(
        [&]()
        {
            out.generate_r1cs_constraints();
            trans.generate_r1cs_constraints();
            for (size_t i = 0; i < other.size(); ++i)
                for (size_t j = 0; j < other[i].size(); ++j)
                    other[i][j].generate_r1cs_constraints();
            gadget->generate_r1cs_constraints();
        },
        1, 1, "Constraint generation", false);
    log_file << elap << '\t';
    log_file.flush();

    trans.generate_r1cs_witness(tree->get_node(trans_idx)->get_digest());
    pb.val(idx) = trans_idx;

    const Node *aux = tree->get_node(trans_idx)->get_f();
    for (size_t i = 0; i < other.size(); ++i, aux = aux->get_f())
        for (size_t j = 0; j < other[i].size(); ++j)
            other[i][j].generate_r1cs_witness(aux->get_c(j)->get_digest());

    // Witness generation
    elap = measure([&]() { gadget->generate_r1cs_witness(); }, 1, 1, "Witness generation", false);
    log_file << elap << '\t';
    log_file.flush();

    // Key generation
    r1cs_ppzksnark_keypair<ppT> keypair;
    elap = measure(
        [&]() { keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system()); }, 1,
        1, "Key generation", false);
    log_file << elap << '\t';
    log_file.flush();

    // Proof generation
    libsnark::r1cs_ppzksnark_proof<ppT> proof;
    elap = measure(
        [&]()
        {
            proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                         pb.auxiliary_input());
        },
        1, 1, "Proof generation", false);
    log_file << elap << '\t';
    log_file.flush();

    // Proof Verification
    bool result;
    elap = measure(
        [&]()
        {
            result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk,
                                                                      pb.primary_input(), proof);
        },
        1, 1, "Proof verification", false);
    log_file << elap << '\n';
    log_file.flush();


    return result;
}

template<size_t first, size_t last, typename GadHash>
void test_range(const char *name)
{
    if constexpr (first < last)
    {
        log_file << first << '\t';
        log_file.flush();

        test_mtree<first, GadHash>();

        test_range<first * 2, last, GadHash>(name);
    }
}

int main()
{
    fs::create_directories("./log");
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                               std::chrono::system_clock::now().time_since_epoch())
                                               .count());
    std::string log_file_name = std::string("./log/mtree_benchmark_") + timestamp +
                                std::string(".log");
    log_file.open(log_file_name);

    std::cout << "Logging to " << log_file_name << "...\n";

    log_file << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    libff::default_ec_pp::init_public_params();

    log_file << "Merkle Tree Benchmark"
             << "\n";
    log_file << "Prime:\t" << FieldT::mod << "\n";
    log_file << "d:\t"
             << "5"
             << "\n";
    log_file << "Minimum Merkle Tree Height:\t" << MIN_HEIGHT << "\n";
    log_file << "Maximum Merkle Tree Height:\t" << MAX_HEIGHT - 1 << "\n\n";
    std::string table_header = std::string("Height\t") + std::string("Tree\t") +
                               std::string("Gadget\t") + std::string("Constraint\t") +
                               std::string("Witness\t") + std::string("Key\t") +
                               std::string("Proof\t") + std::string("Verify\n");


    /**
    To benchmark some permutation gadget over a merkle tree, follow the syntax below:

    log_file << "Permutation_Name\n";
    log_file << table_header;
    test_range<min_tree_height, max_tree_height, PermutationGadget>("Permutation_Name"); 
    lof_file << "\n";

    Different gadgets can require different template arguments to be instantiated. Refer to their 
    documentation for more details. 
    Example:
    GriffinGadget, Poseidon5Gadget and ArionGadget require an instantiation of the respective
    permutation as argument. 
    In particular, an instantiation of ArionGadget follows the syntax below:
        PermutationGadget = ArionGadget<Arion<FieldT, rate, capacity, rounds>>
    Griffin follows the same syntax, while Poseidon require both the number of full rounds and the 
    number of partial rounds.
    **/

    log_file << "Griffin (2:1) r = 12\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, GriffinGadget<Griffin<FieldT, 2, 1, 12>>>("Griffin");
    log_file << "\n";

    log_file << "Griffin (3:1) r = 11\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, GriffinGadget<Griffin<FieldT, 3, 1, 11>>>("Griffin");
    log_file << "\n";

    log_file << "Griffin (7:1) r = 9\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, GriffinGadget<Griffin<FieldT, 7, 1, 9>>>("Griffin");
    log_file << "\n";

    log_file << "Poseidon5 (2:1)\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, Poseidon5Gadget<Poseidon5<FieldT, 2, 1, 4, 56>>>(
        "Poseidon5");
    log_file << "\n";

    log_file << "Poseidon5 (3:1) r_f = 8 r_p = 56\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, Poseidon5Gadget<Poseidon5<FieldT, 3, 1, 4, 56>>>(
        "Poseidon5");
    log_file << "\n";

    log_file << "Poseidon5 (7:1) r_f = 8 r_p = 56\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, Poseidon5Gadget<Poseidon5<FieldT, 7, 1, 4, 56>>>(
        "Poseidon5");
    log_file << "\n";

    log_file << "Arion (2:1) r = 10\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 2, 1, 10>>>("Arion");
    log_file << "\n";

    log_file << "Arion (3:1) r = 8\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 3, 1, 8>>>("Arion");
    log_file << "\n";

    log_file << "Arion (7:1) r = 6\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 7, 1, 6>>>("Arion");
    log_file << "\n";

    log_file << "Aggressive Arion (2:1) r = 7\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 2, 1, 7>>>("Aggressive Arion");
    log_file << "\n";

    log_file << "Aggressive Arion (3:1) r = 6\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 3, 1, 6>>>("Aggressive Arion");
    log_file << "\n";

    log_file << "Aggressive Arion (7:1) r = 4\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionGadget<Arion<FieldT, 7, 1, 4>>>("Aggressive Arion");
    log_file << "\n";

    log_file << "ArionV2 (2:1) r = 6\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 2, 1, 6>>>("ArionV2");
    log_file << "\n";

    log_file << "ArionV2 (3:1) r = 5\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 3, 1, 5>>>("ArionV2");
    log_file << "\n";

    log_file << "ArionV2 (7:1) r = 4\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 7, 1, 4>>>("ArionV2");
    log_file << "\n";

    log_file << "Aggressive ArionV2 (2:1) r = 4\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 2, 1, 4>>>("ArionV2");
    log_file << "\n";

    log_file << "Aggressive ArionV2 (3:1) r = 4\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 3, 1, 4>>>("ArionV2");
    log_file << "\n";

    log_file << "Aggressive ArionV2 (7:1) r = 4\n";
    log_file << table_header;
    test_range<MIN_HEIGHT, MAX_HEIGHT, ArionV2Gadget<ArionV2<FieldT, 7, 1, 4>>>("ArionV2");
    log_file << "\n";

    log_file.close();

    return 0;
}
