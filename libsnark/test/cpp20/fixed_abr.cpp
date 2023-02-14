#include "tree/fixed_abr.hpp"
#include "hash/mimc256.hpp"
#include "hash/mimc512f.hpp"
#include "hash/sha256.hpp"
#include "hash/sha512.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

static bool run_tests()
{
    bool check = true;
    bool all_check = true;

    static constexpr size_t HEIGHT = 4;
    static constexpr auto digest_sha256 =
        BIGHEX(e54f319bda1edc07b45f34a5b6452a2c75bee8332a65ecf5c1803534b9b6e372);

    static constexpr auto digest_sha512 =
        BIGHEX(8eb195cebaf15f4a0c277829505d9b4eedf0d0167183fea9ee74ec93eab6192f37d8857b5d8ba5573300357b92142c906eb9b4ffa6f0297f8c538b81865fef0d);


    std::cout << std::boolalpha;

    std::cout << "ABR SHA256... ";
    check = true;
    {
        std::vector<uint8_t> data(FixedAbr<HEIGHT, Sha256>::INPUT_SIZE);
        FixedAbr<HEIGHT, Sha256> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest_sha256.data(), digest_sha256.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "ABR SHA512... ";
    check = true;
    {
        std::vector<uint8_t> data(FixedAbr<HEIGHT, Sha512>::INPUT_SIZE);
        FixedAbr<HEIGHT, Sha512> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest_sha512.data(), digest_sha512.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "ABRPath SHA256... ";
    check = true;
    {
        std::vector<uint8_t> data(FixedAbr<HEIGHT, Sha512>::INPUT_SIZE);
        FixedAbr<HEIGHT, Sha512> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest_sha512.data(), digest_sha512.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;



    return all_check;
}

int main()
{
    std::cout << "\n==== Testing ABR ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
