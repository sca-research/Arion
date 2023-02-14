#pragma once

#include "util/const_math.hpp"
#include "util/string_utils.hpp"

#include <array>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>

#if __cplusplus >= 202002L
    #include <ranges>
#endif

template<typename Hash>
class MTreeNode
{
public:
    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;

private:
    std::array<uint8_t, Hash::DIGEST_SIZE> digest;
    MTreeNode *f;
    std::array<MTreeNode *, ARITY> c;
    size_t depth;

    template<size_t, typename>
    friend class MTree;

    template<size_t, typename>
    friend class MTreePath;

public:
    MTreeNode() = default;

    MTreeNode(const void *digest, size_t depth) : f{nullptr}, c{}, depth{depth}
    {
        memcpy(this->digest.data(), digest, Hash::DIGEST_SIZE);
    }

    MTreeNode(const std::array<const void *, ARITY> &data, size_t depth) :
        f{nullptr}, c{}, depth{depth}
    {
        std::array<uint8_t, Hash::BLOCK_SIZE> block;

        for (size_t i = 0; i < ARITY; ++i)
            memcpy(block.data() + i * Hash::DIGEST_SIZE, data[i], Hash::DIGEST_SIZE);

        Hash::hash_oneblock(this->digest.data(), block.data());
    }

    const auto &get_digest() const { return digest; }
    const MTreeNode *get_f() const { return f; }
    const MTreeNode *get_c(size_t i) const { return c[i]; }

    friend std::ostream &operator<<(std::ostream &os, const MTreeNode &node)
    {
        for (size_t i = 0; i < node.depth; ++i)
            os << "    ";

        os << "*: " << hexdump(node.digest, false, 64) << '\n';

        for (size_t i = 0; i < ARITY; ++i)
            if (node.c[i] != nullptr)
                os << *node.c[i];

        return os;
    }
};

template<size_t height, typename Hash>
class MTree
{
public:
    using Node = MTreeNode<Hash>;

    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;
    static constexpr size_t LEAVES_N = pow(ARITY, height - 1);
    static constexpr size_t NODES_N = pow_sum(ARITY, (size_t)0, height);
    static constexpr size_t INPUT_SIZE = LEAVES_N * Hash::BLOCK_SIZE;

private:
    std::vector<Node> nodes{};
    Node *root = nullptr;

public:
    MTree() = default;
#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    MTree(const Range &range) :
        MTree(std::ranges::cdata(range),
              std::ranges::size(range) * sizeof(*std::ranges::cdata(range)))
    {}
#endif

    template<typename Iter>
    MTree(const Iter begin, const Iter end) :
        MTree(&*begin, std::distance(begin, end) * sizeof(*begin))
    {}

    MTree(const void *vdata, size_t sz) : nodes(NODES_N), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "MTree: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

#pragma omp parallel for
        // add leaves
        for (size_t i = 0; i < LEAVES_N; ++i)
        {
            std::array<const void *, ARITY> children;

            for (size_t j = 0; j < ARITY; ++j)
                children[j] = data + i * Hash::BLOCK_SIZE + j * Hash::DIGEST_SIZE;
            this->nodes[i] = Node{children, depth};
        }

        // build tree bottom-up
        for (size_t i = 0, last = LEAVES_N, len = LEAVES_N; depth > 0; len += pow(ARITY, depth))
        {
            size_t iters = (len - i) / ARITY;
            --depth;
#pragma omp parallel for
            for (size_t j = 0; j < iters; ++j)
            {
                std::array<const void *, ARITY> children;

                for (size_t k = 0; k < ARITY; ++k)
                    children[k] = this->nodes[i + j * ARITY + k].digest.data();

                this->nodes[last + j] = Node{children, depth};

                for (size_t k = 0; k < ARITY; ++k)
                {
                    this->nodes[last + j].c[k] = &this->nodes[i + j * ARITY + k];
                    this->nodes[i + j * ARITY + k].f = &this->nodes[last + j];
                }
            }
            last += iters;
            i += iters * ARITY;
        }
    }

    const uint8_t *digest() const
    {
        return root->digest.data();
    }

    const Node *get_node(size_t i) const
    {
        return &nodes[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const MTree &tree)
    {
        if (!tree.root)
            return os << "*:";

        return os << *tree.root;
    }
};


template<size_t height, typename Hash>
class MTreePath
{
public:
    using Node = MTreeNode<Hash>;

    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;
    static constexpr size_t NODES_N = ARITY * (height - 1) + 1;
    static constexpr size_t INPUT_N = NODES_N - (height - 1);
    static constexpr size_t INPUT_SIZE = INPUT_N * Hash::DIGEST_SIZE;

private:
    std::vector<Node> nodes{};
    Node *root = nullptr;

public:
    MTreePath() = default;
#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    MTreePath(const Range &range) :
        MTreePath(std::ranges::cdata(range),
                  std::ranges::size(range) * sizeof(*std::ranges::cdata(range)))
    {}
#endif

    template<typename Iter>
    MTreePath(const Iter begin, const Iter end) :
        MTreePath(&*begin, std::distance(begin, end) * sizeof(*begin))
    {}

    MTreePath(const void *vdata, size_t sz) : nodes(NODES_N), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "MTreePath: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

        // bootstrap first node of the path
        this->nodes[0] = {data, depth};

        // build tree bottom-up
        for (size_t i = 1; i < height; ++i)
        {
            std::array<const void *, ARITY> children;

            // insert remaining children
            for (size_t j = 1; j < ARITY; ++j)
                this->nodes[(i - 1) * ARITY + j] = {data += Hash::DIGEST_SIZE, depth};

            // collect children
            for (size_t j = 0; j < ARITY; ++j)
                children[j] = this->nodes[(i - 1) * ARITY + j].digest.data();

            --depth;
            // create node and link children
            this->nodes[i * ARITY] = {children, depth};
            for (size_t j = 0; j < ARITY; ++j)
            {
                this->nodes[i * ARITY].c[j] = &this->nodes[(i - 1) * ARITY + j];
                this->nodes[(i - 1) * ARITY + j].f = &this->nodes[i * ARITY];
            }
        }
    }

    const uint8_t *digest() const
    {
        return root->digest.data();
    }

    const Node *get_node(size_t i) const
    {
        return &nodes[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const MTreePath &tree)
    {
        if (!tree.root)
            return os;

        return os << *tree.root;
    }
};
