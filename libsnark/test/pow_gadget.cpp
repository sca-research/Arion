#include "gadget/pow_gadget.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>


using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

bool test(uint64_t p)
{
    FieldT x = FieldT::random_element();
    FieldT y = x ^ p;


    // Test Gadget
    libsnark::protoboard<FieldT> pb;
    PbVariablePP pb_y{pb, FMT("")};
    PbVariablePP pb_x{pb, FMT("")};

    pb.set_input_sizes(1);

    PowGadget<FieldT> gadget{pb, pb_x, p, pb_y, FMT("gadget")};

    gadget.generate_r1cs_constraints();

    pb.val(pb_x) = x;
    gadget.generate_r1cs_witness();

    bool result = pb.val(pb_y) == y;
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


    std::cout << "Power 0... ";
    std::cout.flush();
    check = test(0);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 1... ";
    std::cout.flush();
    check = test(1);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 2... ";
    std::cout.flush();
    check = test(2);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 3... ";
    std::cout.flush();
    check = test(3);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 4... ";
    std::cout.flush();
    check = test(4);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 5... ";
    std::cout.flush();
    check = test(5);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 7... ";
    std::cout.flush();
    check = test(7);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power 2^20... ";
    std::cout.flush();
    check = test(1 << 20);
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Power random... ";
    std::cout.flush();
    check = test(std::random_device{}() % (1 << 20));
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Power Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
