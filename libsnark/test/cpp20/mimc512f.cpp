#include "hash/mimc512f.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using ppT = libff::default_ec_pp;
using FieldT = libff::Fr<ppT>;
using Hash = Mimc512F<FieldT>;

static bool run_tests()
{
    uint8_t msg[Hash::BLOCK_SIZE]{};
    uint8_t dig[Hash::DIGEST_SIZE]{};
    auto real_dig = BIGHEX(0dd86c5a64a28127d248a1bb324baa0aeb682ed509d58a0bc7f636aa515a45b1149883c6c92e4fe4c09cc3e69c895d6bf368a7e54019c602e979d9b3ddbed9bf);
    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;


    std::cout << "Hashing... ";
    check = true;

    Hash::hash_oneblock(dig, msg);
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;

    std::cout << hexdump(dig) << '\n';
    std::cout << check << '\n';
    all_check &= check;

    return all_check;

   return true;
}

int main()
{
    std::cout << "\n==== Testing MiMC512F ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
