#include "hash/sha256.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

static bool run_tests()
{
    uint8_t msg[Sha256::BLOCK_SIZE]{0x80};
    uint8_t dig[Sha256::DIGEST_SIZE]{};
    auto real_dig = BIGHEX(e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855);
    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;

    std::cout << "Hashing... ";
    check = true;

    Sha256::hash_oneblock(dig, msg);
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;

    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing SHA256 ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
