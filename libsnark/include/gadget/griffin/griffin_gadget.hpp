#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"
#include "gadget/pb_variable_pp.hpp"
#include "hash/griffin.hpp"

template<typename Griffin>
class GriffinGadget : public GadgetPP<typename Griffin::Field>
{
public:
    using Field = typename Griffin::Field;
    using Hash = Griffin;
    using DigVar = FieldVariable<Field>;
    using BlockVar = std::array<DigVar, Hash::RATE>;

    static constexpr size_t RATE = Hash::RATE;
    static constexpr size_t CAPACITY = Hash::CAPACITY;
    static constexpr size_t BRANCH_N = Hash::BRANCH_N;
    static constexpr size_t ROUNDS_N = Hash::ROUNDS_N;
    static constexpr size_t DIGEST_SIZE = Hash::DIGEST_SIZE;
    static constexpr size_t BLOCK_SIZE = Hash::BLOCK_SIZE;
    static constexpr size_t DIGEST_VARS = 1;

    const BlockVar in;
    const DigVar out;

private:
    using super = GadgetPP<typename Griffin::Field>;
    using LC = typename super::LC;

    using super::constrain;
    using super::val;

    static constexpr size_t INTER0_N = 3 * ROUNDS_N;
    static constexpr size_t INTER1_N = 3 * ROUNDS_N;
    static constexpr size_t INTERk_N = 2 * ROUNDS_N;
    static constexpr size_t BRANCH_N4 = BRANCH_N / 4;

    static constexpr auto &alpha = Hash::alpha;
    static constexpr auto &gamma = Hash::gamma;
    static constexpr auto &rc = Hash::round_c;
    static constexpr auto &circ_mat = Hash::circ_mat;

    std::vector<PbVariablePP<Field>> inter[BRANCH_N];

public:
    GriffinGadget(libsnark::protoboard<Field> &pb, const BlockVar &in, const DigVar &out,
                  const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        in{in}, out{out}
    {
        for (size_t i = 0; i < INTER0_N; ++i)
            inter[0].emplace_back(pb, FMT(""));

        for (size_t i = 0; i < INTER1_N; ++i)
            inter[1].emplace_back(pb, FMT(""));

        for (size_t i = 2; i < BRANCH_N; ++i)
            for (size_t j = 0; j < INTERk_N; ++j)
                inter[i].emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {
        size_t i[BRANCH_N]{};
        LC t[BRANCH_N]{};
        LC s[BRANCH_N]{};
        LC l;

        for (size_t j = 0; j < RATE; ++j)
            s[j] = in[j][0];

        // first circulant matrix
        if constexpr (BRANCH_N <= 4)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                for (size_t k = 0; k < BRANCH_N; ++k)
                    t[j] = t[j] + circ_mat[(BRANCH_N - j + k) % BRANCH_N] * s[k];
        }
        else
        {
            static_assert(BRANCH_N % 4 == 0, "BRANCH_N must be a multiple of 4");

            for (size_t j = 0; j < BRANCH_N4; ++j)
                for (size_t k1 = 0; k1 < BRANCH_N4; ++k1)
                    for (size_t k2 = 0, off = 4 * (j != k1); k2 < 4; ++k2)
                        for (size_t k3 = 0; k3 < 4; ++k3)
                            t[4 * j + k2] = t[4 * j + k2] + circ_mat[off + k3] * s[4 * k1 + k3];
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // Base case, y[0] = x[0]^e = x[0]^(1/d) (d = 5)
            i[0] += constrain(inter[0][i[0]], inter[0][i[0] + 2], t[0]);
            i[0] += constrain(inter[0][i[0]], inter[0][i[0]], inter[0][i[0] - 1]);
            i[0] += constrain(inter[0][i[0]], inter[0][i[0]], inter[0][i[0] - 1]);

            // Base case, y[1] = x[1]^d (d = 5)
            i[1] += constrain(t[1], t[1], inter[1][i[1]]);
            i[1] += constrain(inter[1][i[1] - 1], inter[1][i[1] - 1], inter[1][i[1]]);
            i[1] += constrain(inter[1][i[1] - 1], t[1], inter[1][i[1]]);

            // Recursive case y[i] = x[i] * (L(y0,y1,old)^2 + a1*L(y0,y1,old) + a2)
            // <==> y[i] = x[i] * (L(y0,y1,old) * (L(y0,y1,old) + a1) + a2)
            // L(y1, y2, old) = gamma*y1 + y2 + old <==> L(x1, x2, 0) = gamma*x1 + x2
            for (size_t k = 2; k < BRANCH_N; ++k)
            {
                l = gamma * inter[0][i[0] - 1] + inter[1][i[1] - 1];
                if (k != 2)
                    l = l + t[k - 1];
                // y = x*(l^2 + a*l + b) <==> y' - (a*l + b) = l^2 && y = x*y'
                i[k] += constrain(l, l, inter[k][i[k]] - (alpha.first * l + alpha.second));
                i[k] += constrain(t[k], inter[k][i[k] - 1], inter[k][i[k]]);
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // CIRCULANT MATRIX
            if constexpr (BRANCH_N <= 4)
            {
                for (size_t k1 = 0; k1 < BRANCH_N; ++k1)
                    for (size_t k2 = 0; k2 < BRANCH_N; ++k2)
                        t[k1] = t[k1] +
                                circ_mat[(BRANCH_N - k1 + k2) % BRANCH_N] * inter[k2][i[k2] - 1];
            }
            else
            {
                for (size_t k1 = 0; k1 < BRANCH_N4; ++k1)
                    for (size_t k2 = 0; k2 < BRANCH_N4; ++k2)
                        for (size_t k3 = 0, off = 4 * (k1 != k2); k3 < 4; ++k3)
                            for (size_t k4 = 0; k4 < 4; ++k4)
                                t[4 * k1 + k3] = t[4 * k1 + k3] +
                                                 circ_mat[off + k4] *
                                                     inter[4 * k2 + k4][i[4 * k2 + k4] - 1];
            }
            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = t[k] + rc[j * BRANCH_N + k];
        }

        constrain(t[0], 1, out[0]);
    }

    void generate_r1cs_witness()
    {
        size_t i[BRANCH_N]{};
        Field s[BRANCH_N]{};
        Field t[BRANCH_N]{};
        Field l;

        for (size_t j = 0; j < RATE; ++j)
            s[j] = val(in[j][0]);

        // first circulant matrix
        if constexpr (BRANCH_N <= 4)
        {
            for (size_t j = 0; j < BRANCH_N; ++j)
                for (size_t k = 0; k < BRANCH_N; ++k)
                    t[j] = t[j] + circ_mat[(BRANCH_N - j + k) % BRANCH_N] * s[k];
        }
        else
        {
            static_assert(BRANCH_N % 4 == 0, "BRANCH_N must be a multiple of 4");

            for (size_t j = 0; j < BRANCH_N4; ++j)
                for (size_t k1 = 0; k1 < BRANCH_N4; ++k1)
                    for (size_t k2 = 0, off = 4 * (j != k1); k2 < 4; ++k2)
                        for (size_t k3 = 0; k3 < 4; ++k3)
                            t[4 * j + k2] += circ_mat[off + k3] * s[4 * k1 + k3];
        }

        for (size_t j = 0; j < ROUNDS_N; ++j)
        {
            // GTDS
            // y = x^(1/5)
            val(inter[0][i[0] + 2]) = t[0];
            Hash::fifth_inv(val(inter[0][i[0] + 2]));
            val(inter[0][i[0] + 1]) = val(inter[0][i[0] + 2]) * val(inter[0][i[0] + 2]);
            val(inter[0][i[0]]) = val(inter[0][i[0] + 1]) * val(inter[0][i[0] + 1]);
            i[0] += 3;

            // x^5
            val(inter[1][i[1]]) = t[1] * t[1];
            ++i[1];
            val(inter[1][i[1]]) = val(inter[1][i[1] - 1]) * val(inter[1][i[1] - 1]);
            ++i[1];
            val(inter[1][i[1]]) = val(inter[1][i[1] - 1]) * t[1];
            ++i[1];

            // L(x1, x2, 0) = gamma*x1 + x2
            l = gamma * val(inter[0][i[0] - 1]) + val(inter[1][i[1] - 1]);

            // Recursive case y[i] = x[i] * (L(y0,y1,old)^2 + a1*L(y0,y1,old) + a2)
            for (size_t k = 2; k < BRANCH_N; ++k)
            {
                l = gamma;
                l *= val(inter[0][i[0] - 1]);
                l += val(inter[1][i[1] - 1]);
                if (k != 2)
                    l += t[k - 1];
                val(inter[k][i[k]]) = l;
                val(inter[k][i[k]]) += alpha.first;
                val(inter[k][i[k]]) *= l;
                val(inter[k][i[k]]) += alpha.second;
                ++i[k];
                val(inter[k][i[k]]) = val(inter[k][i[k] - 1]);
                val(inter[k][i[k]]) *= t[k];
                ++i[k];
            }

            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] = 0;

            // CIRCULANT MATRIX
            if constexpr (BRANCH_N <= 4)
            {
                for (size_t k1 = 0; k1 < BRANCH_N; ++k1)
                    for (size_t k2 = 0; k2 < BRANCH_N; ++k2)
                        t[k1] += circ_mat[(BRANCH_N - k1 + k2) % BRANCH_N] *
                                 val(inter[k2][i[k2] - 1]);
            }
            else
            {
                for (size_t k1 = 0; k1 < BRANCH_N4; ++k1)
                    for (size_t k2 = 0; k2 < BRANCH_N4; ++k2)
                        for (size_t k3 = 0, off = 4 * (k1 != k2); k3 < 4; ++k3)
                            for (size_t k4 = 0; k4 < 4; ++k4)
                                t[4 * k1 + k3] += circ_mat[off + k4] *
                                                  val(inter[4 * k2 + k4][i[4 * k2 + k4] - 1]);
            }

            // ADD CONSTANTS
            for (size_t k = 0; k < BRANCH_N; ++k)
                t[k] += rc[j * BRANCH_N + k];
        }

        val(out[0]) = t[0];
    }
};
