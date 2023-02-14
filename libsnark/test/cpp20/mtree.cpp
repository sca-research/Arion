#include "tree/mtree.hpp"
#include "hash/sha256.hpp"
#include "hash/sha512.hpp"
#include "hash/arion.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using FieldT = libff::Fq<libff::default_ec_pp>;

static bool run_tests()
{
    bool check = true;
    bool all_check = true;

    static constexpr size_t HEIGHT = 4;

    auto digest256 = BIGHEX(26b0052694fc42fdff93e6fb5a71d38c3dd7dc5b6ad710eb048c660233137fab);
    auto digest512 =
        BIGHEX(6e3d539e81fcba88a5a6875590df1f6ec06e67b1656504bd60f33953b81b080637d790e65789a4e03bfa4d457cb820f5153a0299c74798775d4295e9b0955517);

    std::cout << std::boolalpha;


    std::cout << "Full Tree SHA256... ";
    check = true;
    {
        std::vector<uint8_t> data(MTree<HEIGHT, Sha256>::INPUT_SIZE);
        MTree<HEIGHT, Sha256> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest256.data(), digest256.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Full Tree SHA512... ";
    check = true;
    {
        std::vector<uint8_t> data(MTree<HEIGHT, Sha512>::INPUT_SIZE);
        MTree<HEIGHT, Sha512> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest512.data(), digest512.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Full Tree Arion... ";
    check = true;
    {
        std::vector<uint8_t> data(MTree<HEIGHT, Arion<FieldT, 2, 1>>::INPUT_SIZE);
        MTree<HEIGHT, Arion<FieldT, 2, 1>> tree(data.begin(), data.end());

        std::cout << '\n' << tree << '\n';
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Tree Path SHA256... ";
    check = true;
    {
        std::vector<uint8_t> data(MTreePath<HEIGHT, Sha256>::INPUT_SIZE);
        uint8_t block[Sha256::BLOCK_SIZE]{};

        Sha256::hash_oneblock(data.data(), block);
        Sha256::hash_oneblock(data.data() + Sha256::DIGEST_SIZE, block);

        for (size_t i = Sha256::BLOCK_SIZE; i < data.size(); i += Sha256::DIGEST_SIZE)
        {
            memcpy(block, data.data() + i - Sha256::DIGEST_SIZE, Sha256::DIGEST_SIZE);
            memcpy(block + Sha256::DIGEST_SIZE, data.data() + i - Sha256::DIGEST_SIZE,
                   Sha256::DIGEST_SIZE);

            Sha256::hash_oneblock(data.data() + i, block);
        }

        MTreePath<HEIGHT, Sha256> tree(data.begin(), data.end());

        std::cout << '\n' << tree << '\n';
        check = memcmp(tree.digest(), digest256.data(), digest256.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Tree Path SHA512... ";
    check = true;
    {
        std::vector<uint8_t> data(MTreePath<HEIGHT, Sha512>::INPUT_SIZE);
        uint8_t block[Sha512::BLOCK_SIZE]{};

        Sha512::hash_oneblock(data.data(), block);
        Sha512::hash_oneblock(data.data() + Sha512::DIGEST_SIZE, block);

        for (size_t i = 2 * Sha512::DIGEST_SIZE; i < data.size(); i += Sha512::DIGEST_SIZE)
        {
            memcpy(block, data.data() + i - Sha512::DIGEST_SIZE, Sha512::DIGEST_SIZE);
            memcpy(block + Sha512::DIGEST_SIZE, data.data() + i - Sha512::DIGEST_SIZE,
                   Sha512::DIGEST_SIZE);

            Sha512::hash_oneblock(data.data() + i, block);
        }

        MTreePath<HEIGHT, Sha512> tree(data.begin(), data.end());

        check = memcmp(tree.digest(), digest512.data(), digest512.size()) == 0;
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Tree path Arion... ";
    check = true;
    {
        std::vector<uint8_t> data(MTreePath<HEIGHT, Arion<FieldT, 2, 1>>::INPUT_SIZE);
        MTreePath<HEIGHT, Arion<FieldT, 2, 1>> tree(data.begin(), data.end());

        std::cout << '\n' << tree << '\n';
    }
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Merkle Tree ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
