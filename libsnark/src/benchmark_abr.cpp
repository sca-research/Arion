#include "gadget/abr_gadget.hpp"
#include "tree/fixed_abr.hpp"
#include "util/measure.hpp"

#include "gadget/griffin/griffin_gadget.hpp"
#include "gadget/mimc256/mimc256_gadget.hpp"
#include "gadget/mimc512f/mimc512f_gadget.hpp"
#include "gadget/mimc512f2k/mimc512f2k_gadget.hpp"
#include "gadget/sha256/sha256_gadget.hpp"
#include "gadget/sha512/sha512_gadget.hpp"
#include "gadget/poseidon5/poseidon5_gadget.hpp"
#include "gadget/arion/arion_gadget.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>
#include <omp.h>
#include <chrono>
#include <string> 
#include <filesystem>
namespace fs = std::filesystem;

static constexpr size_t TRANS_IDX = 0;
static constexpr size_t MIN_TREE_HEIGHT = 4;
static constexpr size_t MAX_TREE_HEIGHT = 31;

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

using GadGriffin = griffin_two_to_one_hash_gadget<FieldT>;
using GadSha256 = libsnark::sha256_two_to_one_hash_gadget<FieldT>;
using GadSha512 = libsnark::sha512::sha512_two_to_one_hash_gadget<FieldT>;
using GadMimc256 = mimc256_two_to_one_hash_gadget<FieldT>;
using GadMimc512F = mimc512f_two_to_one_hash_gadget<FieldT>;
using GadMimc512F2K = mimc512f2k_two_to_one_hash_gadget<FieldT>;
using GadPoseidon5 = poseidon5_two_to_one_hash_gadget<FieldT>;
using GadArion = arion_two_to_one_hash_gadget<FieldT>;

std::ofstream log_file;

template<size_t tree_height, typename GadHash>
bool test_tRee()
{
    using Hash = typename GadHash::Hash;
    using GadTree = ABR_Gadget<FieldT, GadHash>;
    using DigVar = typename GadTree::DigVar;
    using FixTree = FixedAbr<tree_height, Hash>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    double elap = 0;

    // Build tree
    libff::bit_vector trans_bv;
    libff::bit_vector other_bv;
    std::vector<libff::bit_vector> middle_bv;
    std::vector<libff::bit_vector> otherx_bv;
    libff::bit_vector out_bv;

    elap = measure(
        [&]()
        {
            std::vector<uint8_t> data(FixTree::INPUT_SIZE);
            std::generate(data.begin(), data.end(), std::ref(rng));
            FixTree tree{data.begin(), data.end()};

            // Extract our transaction
            trans_bv.resize(DIGEST_VARS);
            unpack_bits(trans_bv, tree.get_node(TRANS_IDX)->get_digest());

            // Extract other transaction node
            other_bv.resize(DIGEST_VARS);
            unpack_bits(other_bv, tree.get_node(TRANS_IDX + 1)->get_digest());

            // Extract middle nodes
            middle_bv.resize(tree_height - 2);
            std::fill(middle_bv.begin(), middle_bv.end(), libff::bit_vector(DIGEST_VARS));
            for (size_t i = 0, j = FixTree::LEAVES_N; i < middle_bv.size();
                 ++i, j += 1ULL << (tree_height - 2 - i))
                unpack_bits(middle_bv[i], tree.get_node(j)->get_digest());

            // Extract otherx nodes
            otherx_bv.resize(tree_height - 2);
            std::fill(otherx_bv.begin(), otherx_bv.end(), libff::bit_vector(DIGEST_VARS));
            for (size_t i = 0, j = FixTree::INPUT_N + 1; i < otherx_bv.size();
                 ++i, j += 1ULL << (tree_height - 1 - i))
                unpack_bits(otherx_bv[i], tree.get_node(j)->get_digest());

            // Extract output node
            out_bv.resize(DIGEST_VARS);
            unpack_bits(out_bv, tree.digest());
        },
        1, 1, "Tree Generation", false);
    log_file << elap << '\t';
    log_file.flush();

    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    std::vector<DigVar> out;
    std::vector<DigVar> trans;
    std::vector<DigVar> other;
    std::vector<DigVar> middle;
    std::vector<DigVar> otherx;
    std::vector<GadTree> gadget;

    elap = measure(
        [&]()
        {
            out.emplace_back(pb, DIGEST_VARS, FMT("out"));
            trans.emplace_back(pb, DIGEST_VARS, FMT("trans"));
            other.emplace_back(pb, DIGEST_VARS, FMT("other"));

            for (size_t i = 0; i < middle_bv.size(); ++i)
                middle.emplace_back(pb, DIGEST_VARS, FMT("middle_%llu", i));

            for (size_t i = 0; i < otherx_bv.size(); ++i)
                otherx.emplace_back(pb, DIGEST_VARS, FMT("otherx_%llu", i));

            pb.set_input_sizes(DIGEST_VARS);
            gadget.emplace_back(pb, out[0], trans[0], other[0], middle, otherx, TRANS_IDX,
                                tree_height, FMT("merkle_tree"));
        },
        1, 1, "Gadget construction", false);
    log_file << elap << '\t';
    log_file.flush();

    elap = measure(
        [&]()
        {
            out[0].generate_r1cs_constraints();
            trans[0].generate_r1cs_constraints();
            other[0].generate_r1cs_constraints();
            for (auto &&x : middle)
                x.generate_r1cs_constraints();
            for (auto &&x : otherx)
                x.generate_r1cs_constraints();
            gadget[0].generate_r1cs_constraints();
        },
        1, 1, "Constraint generation", false);
    log_file << elap << '\t';
    log_file.flush();

    elap = measure(
        [&]()
        {
            out[0].generate_r1cs_witness(out_bv);
            trans[0].generate_r1cs_witness(trans_bv);
            other[0].generate_r1cs_witness(other_bv);
            for (size_t i = 0; i < middle.size(); ++i)
                middle[i].generate_r1cs_witness(middle_bv[i]);
            for (size_t i = 0; i < otherx.size(); ++i)
                otherx[i].generate_r1cs_witness(otherx_bv[i]);
            gadget[0].generate_r1cs_witness();
        },
        1, 1, "Witness generation", false);
    log_file << elap << '\t';
    log_file.flush();

    bool result;
    std::vector<libsnark::r1cs_ppzksnark_keypair<ppT>> keypair;
    libsnark::r1cs_ppzksnark_proof<ppT> proof;

    elap = measure(
        [&]() {
            keypair.emplace_back(
                libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system()));
        },
        1, 1, "Key generation", false);
    log_file << elap << '\t';
    log_file.flush();

    elap = measure(
        [&]()
        {
            proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair[0].pk, pb.primary_input(),
                                                         pb.auxiliary_input());
        },
        1, 1, "Proof generation", false);
    log_file << elap << '\t';
    log_file.flush();

    elap = measure(
        [&]()
        {
            result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair[0].vk,
                                                                      pb.primary_input(), proof);
        },
        1, 1, "Proof verification", false);
    log_file << elap << '\n';
    log_file.flush();


    return result;
}

template<size_t tree_height, typename GadHash>
bool test_ptRee()
{
    using Hash = typename GadHash::Hash;
    using GadTree = ABR_Gadget<FieldT, GadHash>;
    using DigVar = typename GadTree::DigVar;
    using FixTree = FixedAbr<tree_height, Hash>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};
    double elap = 0;

    // Build tree
    FixTree tree;
    elap = measure(
        [&]()
        {
            std::vector<uint8_t> data(FixTree::INPUT_SIZE);
            std::generate(data.begin(), data.end(), std::ref(rng));
            tree = FixTree{data.begin(), data.end()};
        },
        1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // There is no need to extract anything as FieldVariable does that for us

    // Build Gadget
    libsnark::protoboard<FieldT> pb;
    std::vector<DigVar> out;
    std::vector<DigVar> trans;
    std::vector<DigVar> other;
    std::vector<DigVar> middle;
    std::vector<DigVar> otherx;
    std::vector<GadTree> gadget;

    elap = measure(
        [&]()
        {
            out.emplace_back(pb, DIGEST_VARS, FMT("out"));
            trans.emplace_back(pb, DIGEST_VARS, FMT("trans"));
            other.emplace_back(pb, DIGEST_VARS, FMT("other"));
            for (size_t i = 0; i < tree_height - 2; ++i)
                middle.emplace_back(pb, DIGEST_VARS, FMT("middle_%llu", i));
            for (size_t i = 0; i < tree_height - 2; ++i)
                otherx.emplace_back(pb, DIGEST_VARS, FMT("otherx_%llu", i));


            pb.set_input_sizes(1);
            gadget.emplace_back(pb, out[0], trans[0], other[0], middle, otherx, TRANS_IDX,
                                tree_height, FMT("merkle_tree"));
        },
        1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // Generate constraints
    elap = measure([&]() { gadget[0].generate_r1cs_constraints(); }, 1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // Generate witnesses
    elap = measure(
        [&]()
        {
            out[0].generate_r1cs_witness(tree.digest());
            trans[0].generate_r1cs_witness(tree.get_node(TRANS_IDX)->get_digest());
            other[0].generate_r1cs_witness(tree.get_node(TRANS_IDX + 1)->get_digest());

            for (size_t i = 0, j = FixTree::LEAVES_N; i < middle.size();
                 ++i, j += 1ULL << (tree_height - 2 - i))
                middle[i].generate_r1cs_witness(tree.get_node(j)->get_digest());

            for (size_t i = 0, j = FixTree::INPUT_N + 1; i < otherx.size();
                 ++i, j += 1ULL << (tree_height - 1 - i))
                otherx[i].generate_r1cs_witness(tree.get_node(j)->get_digest());

            gadget[0].generate_r1cs_witness();
        },
        1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // Generate key
    bool result;
    std::vector<libsnark::r1cs_ppzksnark_keypair<ppT>> keypair;
    elap = measure(
        [&]() {
            keypair.emplace_back(
                libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system()));
        },
        1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // Generate proof
    libsnark::r1cs_ppzksnark_proof<ppT> proof;
    elap = measure(
        [&]()
        {
            proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair[0].pk, pb.primary_input(),
                                                         pb.auxiliary_input());
        },
        1, 1, "", false);
    log_file << elap << '\t';
    log_file.flush();

    // Verify proof
    elap = measure(
        [&]()
        {
            result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair[0].vk,
                                                                      pb.primary_input(), proof);
        },
        1, 1, "", false);
    log_file << elap << '\n';
    log_file.flush();


    return result;
}

template<size_t first, size_t last, typename GadHash>
void test_tRee_from(const char *name)
{
    if constexpr (first < last)
    {
        log_file << first << '\t';
        log_file.flush();

        test_tRee<first, GadHash>();

        test_tRee_from<first + 1, last, GadHash>(name);
    }
}

template<size_t first, size_t last, typename GadHash>
void test_ptRee_from(const char *name)
{
    if constexpr (first < last)
    {
        log_file << first << '\t';
        log_file.flush();

        test_ptRee<first, GadHash>();

        test_ptRee_from<first + 1, last, GadHash>(name);
    }
}

int main()
{
    fs::create_directories("./log");
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    std::string log_file_name = std::string("./log/abr_benchmark_") + timestamp + std::string(".log");
    log_file.open(log_file_name);

    log_file << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    log_file << "ABR Benchmark" << "\n";
    log_file << "Prime:\t" << FieldT::mod << "\n";
    log_file << "Minimum Merkle Tree Height:\t" << MIN_TREE_HEIGHT << "\n";
    log_file << "Maximum Merkle Tree Height:\t" << MAX_TREE_HEIGHT - 1 << "\n\n";
    std::string table_header = std::string("Height\t")      +
                               std::string("Tree\t")        +
                               std::string("Gadget\t")      +
                               std::string("Constraint\t")  +
                               std::string("Witness\t")     +
                               std::string("Key\t")         +
                               std::string("Proof\t")       +
                               std::string("Verify\n");

    ppT::init_public_params();
    /*
    log_file << "SHA256\n";
    log_file << table_header;
    test_tRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadSha256>("SHA256");
    log_file << "\n";
    */

    log_file << "MiMC256\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadMimc256>("MiMC256");
    log_file << "\n";*

    log_file << "MiMC512-SK\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadMimc512F>("MiMC512-SK");
    log_file << "\n";*

    log_file << "MiMC512-DK\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadMimc512F2K>("MiMC512-DK");
    log_file << "\n";*

/*
    log_file << "Griffin\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadGriffin>("Griffin");
    log_file << "\n";

    log_file << "Poseidon5\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadPoseidon5>("Poseidon5");
    log_file << "\n";*/

/*
    log_file << "Arion\n";
    log_file << table_header;
    test_ptRee_from<MIN_TREE_HEIGHT, MAX_TREE_HEIGHT, GadArion>("Arion5");
    log_file << "\n";*/

    return 0;
}
