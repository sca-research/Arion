#pragma once

#include "util/string_utils.hpp"

#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>

#if __cplusplus >= 202002L
    #include <ranges>
#endif

template<typename Hash>
class FixedMTreeNode
{
private:
    std::array<uint8_t, Hash::DIGEST_SIZE> digest;
    FixedMTreeNode *f = nullptr;
    FixedMTreeNode *l = nullptr;
    FixedMTreeNode *r = nullptr;
    size_t depth = 0;

    template<size_t sz, typename>
    friend class FixedMTree;

    template<size_t, typename>
    friend class FixedMTreePath;

public:
    FixedMTreeNode() = default;

    FixedMTreeNode(const uint8_t *digest, size_t depth) : depth{depth}
    {
        memcpy(this->digest.data(), digest, Hash::DIGEST_SIZE);
    }

    FixedMTreeNode(const uint8_t *left, const uint8_t *right, size_t depth) : depth{depth}
    {
        uint8_t block[Hash::BLOCK_SIZE]{};

        memcpy(block, left, Hash::DIGEST_SIZE);
        memcpy(block + Hash::DIGEST_SIZE, right, Hash::DIGEST_SIZE);

        Hash::hash_oneblock(this->digest.data(), block);
    }

    const auto &get_digest() const { return digest; }

    friend std::ostream &operator<<(std::ostream &os, const FixedMTreeNode &node)
    {
        for (size_t i = 0; i < node.depth; ++i)
            os << "    ";

        os << "*: " << hexdump(node.get_digest(), false, 64) << '\n';

        if (node.l)
            os << *node.l;
        if (node.r)
            os << *node.r;

        return os;
    }
};


template<size_t height, typename Hash>
class FixedMTree
{
private:
    using Node = FixedMTreeNode<Hash>;

    static constexpr size_t LEAVES_N = 1ULL << (height - 1);

    std::vector<Node> nodes{};
    Node *root = nullptr;

public:
    static constexpr size_t INPUT_SIZE = LEAVES_N * Hash::BLOCK_SIZE;

    FixedMTree() = default;
#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    FixedMTree(const Range &range) :
        FixedMTree(std::ranges::cdata(range),
                   std::ranges::size(range) * sizeof(*std::ranges::cdata(range)))
    {}
#endif

    template<typename Iter>
    FixedMTree(const Iter begin, const Iter end) :
        FixedMTree(&*begin, std::distance(begin, end) * sizeof(*begin))
    {}

    FixedMTree(const void *vdata, size_t sz) : nodes((1ULL << height) - 1), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "FixedMTree: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

// serial code
#if 0
        // add leaves
            for (size_t i = 0; i < LEAVES_N; ++i)
                this->nodes[i] = {data + Hash::BLOCK_SIZE * i,
                                  data + Hash::BLOCK_SIZE * i + Hash::DIGEST_SIZE, depth};

            // build tree bottom-up
            for (size_t i = 0, last = LEAVES_N, len = LEAVES_N; depth > 0; len += 1ULL << depth)
            {
                --depth;
                while (i < len)
                {
                    this->nodes[last] = {this->nodes[i].get_digest().data(), this->nodes[i + 1].get_digest().data(), depth};

                    this->nodes[last].l = &this->nodes[i];
                    this->nodes[last].r = &this->nodes[i + 1];
                    this->nodes[i].f = &this->nodes[last];
                    this->nodes[i + 1].f = &this->nodes[last];
                    i += 2;
                    ++last;
                }
            }
#else // parallel code
    #pragma omp parallel for
        // add leaves
        for (size_t i = 0; i < LEAVES_N; ++i)
            this->nodes[i] = {data + Hash::BLOCK_SIZE * i,
                              data + Hash::BLOCK_SIZE * i + Hash::DIGEST_SIZE, depth};

        // build tree bottom-up
        for (size_t i = 0, last = LEAVES_N, len = LEAVES_N; depth > 0; len += 1ULL << depth)
        {
            size_t iters = (len - i) >> 1;
            --depth;
    #pragma omp parallel for
            for (size_t j = 0; j < iters; ++j)
            {
                size_t k = i + j * 2;
                size_t l = last + j;
                this->nodes[l] = {this->nodes[k].get_digest().data(),
                                  this->nodes[k + 1].get_digest().data(), depth};

                this->nodes[l].l = &this->nodes[k];
                this->nodes[l].r = &this->nodes[k + 1];
                this->nodes[k].f = &this->nodes[l];
                this->nodes[k + 1].f = &this->nodes[l];
            }
            last += iters;
            i += iters * 2;
        }
#endif
    }

    const auto &digest() const
    {
        return root->get_digest();
    }

    const Node *get_node(size_t i) const
    {
        return &nodes[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const FixedMTree &tree)
    {
        if (!tree.root)
            return os;

        return os << *tree.root;
    }
};

template<size_t height, typename Hash>
class FixedMTreePath
{
private:
    using Node = FixedMTreeNode<Hash>;

    static constexpr size_t NODES_N = 2 * height - 1;

    std::vector<Node> nodes{};
    Node *root = nullptr;

public:
    static constexpr size_t INPUT_SIZE = height * Hash::DIGEST_SIZE;

    FixedMTreePath() = default;
#if __cplusplus >= 202002L
    template<std::ranges::range Range>
    FixedMTreePath(const Range &range) :
        FixedMTreePath(std::ranges::cdata(range),
                       std::ranges::size(range) * sizeof(*std::ranges::cdata(range)))
    {}
#endif

    template<typename Iter>
    FixedMTreePath(const Iter begin, const Iter end) :
        FixedMTreePath(&*begin, std::distance(begin, end) * sizeof(*begin))
    {}

    FixedMTreePath(const void *vdata, size_t sz) : nodes(NODES_N), root{&nodes.back()}
    {
        if (sz != INPUT_SIZE)
        {
            std::cerr << "FixedMTreePath: Bad size of input data\n";
            return;
        }

        const uint8_t *data = (const uint8_t *)vdata;
        size_t depth = height - 1;

        // bootstrap first ndoe
        this->nodes[0] = {data, depth};

        // build tree bottom-up
        for (size_t i = 2; i < NODES_N; i += 2)
        {
            // add other children
            this->nodes[i - 1] = {data += Hash::DIGEST_SIZE, depth};

            // build parent
            --depth;
            this->nodes[i] = {this->nodes[i - 2].get_digest().data(),
                              this->nodes[i - 1].get_digest().data(), depth};
            this->nodes[i].l = &this->nodes[i - 2];
            this->nodes[i].r = &this->nodes[i - 1];
            this->nodes[i - 2].f = &this->nodes[i];
            this->nodes[i - 1].f = &this->nodes[i];
        }
    }

    const auto &digest() const
    {
        return root->get_digest();
    }

    const Node *get_node(size_t i) const
    {
        return &nodes[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const FixedMTreePath &tree)
    {
        if (!tree.root)
            return os;

        return os << *tree.root;
    }
};
