#include "hash/sha512.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

static bool run_tests()
{
    uint8_t msg[Sha512::BLOCK_SIZE]{0x80};
    uint8_t dig[Sha512::DIGEST_SIZE]{};
    auto real_dig = BIGHEX(cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e);
    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;

    std::cout << "Hashing... ";
    check = true;

    Sha512::hash_oneblock(dig, msg);
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;

    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing SHA512 ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
