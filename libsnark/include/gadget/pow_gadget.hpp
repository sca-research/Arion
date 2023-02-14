#pragma once

#include "gadget/field_variable.hpp"
#include "gadget/gadget_pp.hpp"
#include "gadget/pb_variable_pp.hpp"

template<typename FieldT>
class PowGadget : public GadgetPP<FieldT>
{
public:
    using super = GadgetPP<FieldT>;
    using PbVar = libsnark::pb_variable<FieldT>;
    using PbLC = PbLCPP<FieldT>;

    using super::constrain;
    using super::val;

    const PbLC x;
    const uint64_t y;
    const PbVar out;

private:
    std::vector<PbVariablePP<FieldT>> inter;
    size_t inter_n;
    uint64_t y_inv;

public:
    PowGadget(libsnark::protoboard<FieldT> &pb, const PbLC &x, uint64_t y, const PbVar &out,
              const std::string &annotation_prefix) :
        super{pb, annotation_prefix},
        x{x}, y{y}, out{out}, inter_n{0}, y_inv{0}
    {
        auto p{y};

        // we need log_2(p) - 2 intermediates (we already have x^1, and we put the result in out)
        while (p > 1)
        {
            inter_n += p % 2;
            y_inv |= p % 2;
            y_inv *= 2;
            p /= 2;
            ++inter_n;
        }
        y_inv += p % 2;
        inter_n -= !!inter_n; // avoid wrap-around

        for (size_t i = 0; i < inter_n; ++i)
            inter.emplace_back(pb, FMT(""));
    }

    void generate_r1cs_constraints()
    {

        switch (y)
        {
        case 0: constrain(1, 1, out); return;
        case 1: constrain(x, 1, out); return;
        case 2: constrain(x, x, out); return;
        case 3:
            constrain(x, x, inter[0]);
            constrain(inter[0], x, out);
            return;
        default:
        {
            size_t i = 1;
            auto p{y};

            constrain(x, x, inter[0]);

            if (p % 2)
            {
                i += constrain(inter[inter_n - i], x, out);
                i += constrain(inter[inter_n - i], inter[inter_n - i], inter[inter_n - i + 1]);
            }
            else
            {
                i += constrain(inter[inter_n - i], inter[inter_n - i], out);
            }
            p /= 2;

            while (i < inter_n)
            {
                if (p % 2)
                    i += constrain(inter[inter_n - i], x, inter[inter_n - i + 1]);
                i += constrain(inter[inter_n - i], inter[inter_n - i], inter[inter_n - i + 1]);
                p /= 2;
            }

            break;
        }
        }
    }

    void generate_r1cs_witness()
    {
        switch (y)
        {
        case 0: val(out) = 1; break;
        case 1: val(out) = val(x); break;
        case 2: val(out) = val(x) * val(x); break;
        default:
        {
            auto p = y_inv;

            val(inter[0]) = val(x) * val(x);

            p /= 2;

            for (size_t i = 1; i < inter_n; ++i)
            {
                if (p % 2)
                {
                    val(inter[i]) = val(inter[i - 1]) * val(x);
                    ++i;
                }
                val(inter[i]) = val(inter[i - 1]) * val(inter[i - 1]);
                p /= 2;
            }

            if (p % 2)
                val(out) = val(inter[inter_n - 1]) * val(x);
            else
                val(out) = val(inter[inter_n - 1]) * val(inter[inter_n - 1]);

            break;
        }
        }

        val(out).as_bigint().print_hex();
    }
};
