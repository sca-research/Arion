#include "hash/sha256.hpp"
#include "util/string_utils.hpp"
#include <iostream>

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Syntax: " << argv[0] << " <block>\n";
        exit(EXIT_FAILURE);
    }

    std::vector<uint8_t> data(Sha256::BLOCK_SIZE);
    std::vector<uint8_t> digest(Sha256::DIGEST_SIZE);

    hexstring_to_range(argv[1], data);

    Sha256::hash_oneblock(digest.data(), data.data());

    std::cout << hexdump(digest) << '\n';

    return 0;
}
