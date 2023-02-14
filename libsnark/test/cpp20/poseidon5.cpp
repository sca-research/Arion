#include "hash/poseidon5.hpp"
#include "util/measure.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using ppT = libff::default_ec_pp;
using FieldT = libff::Fr<ppT>;
using Hash = Poseidon5<FieldT, 2, 1>;

static bool run_tests()
{
    auto msg =
        BIGHEX(00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
    uint8_t dig[Hash::DIGEST_SIZE]{};
    auto real_dig = BIGHEX(115cc0f5e7d690413df64c6b9662e9cf2a3617f2743245519e19607a4417189a);

    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;

    std::cout << "Hashing... ";
    check = true;

    Hash::hash_oneblock(dig, msg.data());
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;

    std::cout << hexdump(dig) << '\n';
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Hash ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
