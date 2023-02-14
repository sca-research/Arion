#include "util/string_utils.hpp"
#include <assert.h>
#include <iostream>

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Syntax: " << argv[0] << " <block_1> <block_2>\n";
        exit(EXIT_FAILURE);
    }

    std::string_view x_str{argv[1]};
    std::string_view y_str{argv[2]};

    if (x_str.size() != y_str.size() || (x_str.size() & 1) || (y_str.size() & 1))
    {
        std::cerr << "Invalid block sizes!\n";
        exit(EXIT_FAILURE);
    }

    size_t sz = x_str.size() / 2;
    std::vector<uint8_t> x(sz);
    std::vector<uint8_t> y(sz);

    hexstring_to_range(x_str, x);
    hexstring_to_range(y_str, y);

    for (size_t i = 0; i < sz; ++i)
        x[i] ^= y[i];

    std::cout << hexdump(x) << '\n';

    return 0;
}
