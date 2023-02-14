/** @file
 *****************************************************************************

 Implementation of interfaces for gadgets for the SHA512 message schedule and round function.

 See sha512_components.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SHA512_COMPONENTS_TCC_
#define SHA512_COMPONENTS_TCC_

namespace libsnark::sha512 {

const unsigned long long SHA512_K[80] =  {
   0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
   0x3956c25bf348b538, 0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
   0xd807aa98a3030242, 0x12835b0145706fbe, 0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
   0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235, 0xc19bf174cf692694,
   0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
   0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
   0x983e5152ee66dfab, 0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4,
   0xc6e00bf33da88fc2, 0xd5a79147930aa725, 0x06ca6351e003826f, 0x142929670a0e6e70,
   0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
   0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
   0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30,
   0xd192e819d6ef5218, 0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8,
   0x19a4c116b8d2d0c8, 0x1e376c085141ab53, 0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
   0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
   0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
   0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b,
   0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
   0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
   0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
   0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817,

};

const unsigned long long SHA512_H[8] = {
        0x6a09e667f3bcc908,
        0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b,
        0xa54ff53a5f1d36f1,
        0x510e527fade682d1,
        0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b,
        0x5be0cd19137e2179,
};

template<typename FieldT>
pb_linear_combination_array<FieldT> SHA512_default_IV(protoboard<FieldT> &pb)
{
    pb_linear_combination_array<FieldT> result;
    result.reserve(SHA512_digest_size);

    for (size_t i = 0; i < SHA512_digest_size; ++i)
    {
        int iv_val = (SHA512_H[i / 64] >> (63-(i % 64))) & 1;

        pb_linear_combination<FieldT> iv_element;
        iv_element.assign(pb, iv_val * ONE);
        iv_element.evaluate(pb);

        result.emplace_back(iv_element);
    }

    return result;
}

template<typename FieldT>
sha512_message_schedule_gadget<FieldT>::sha512_message_schedule_gadget(protoboard<FieldT> &pb,
                                                                       const pb_variable_array<FieldT> &M,
                                                                       const pb_variable_array<FieldT> &packed_W,
                                                                       const std::string &annotation_prefix) :
    gadget<FieldT>(pb, annotation_prefix),
    M(M),
    packed_W(packed_W)
{
    W_bits.resize(80);

    pack_W.resize(16);
    for (size_t i = 0; i < 16; ++i)
    {
        W_bits[i] = pb_variable_array<FieldT>(M.rbegin() + (15-i) * 64, M.rbegin() + (16-i) * 64);
        pack_W[i].reset(new packing_gadget<FieldT>(pb, W_bits[i], packed_W[i], FMT(this->annotation_prefix, " pack_W_%zu", i)));
    }

    /* NB: some of those will be un-allocated */
    sigma0.resize(80);
    sigma1.resize(80);
    compute_sigma0.resize(80);
    compute_sigma1.resize(80);
    unreduced_W.resize(80);
    mod_reduce_W.resize(80);

    for (size_t i = 16; i < 80; ++i)
    {
        /* allocate result variables for sigma0/sigma1 invocations */
        sigma0[i].allocate(pb, FMT(this->annotation_prefix, " sigma0_%zu", i));
        sigma1[i].allocate(pb, FMT(this->annotation_prefix, " sigma1_%zu", i));

        /* compute sigma0/sigma1 */
        compute_sigma0[i].reset(new small_sigma_gadget<FieldT>(pb, W_bits[i-15], sigma0[i], 1, 8, 7, FMT(this->annotation_prefix, " compute_sigma0_%zu", i)));
        compute_sigma1[i].reset(new small_sigma_gadget<FieldT>(pb, W_bits[i-2], sigma1[i], 19, 61, 6, FMT(this->annotation_prefix, " compute_sigma1_%zu", i)));

        /* unreduced_W = sigma0(W_{i-15}) + sigma1(W_{i-2}) + W_{i-7} + W_{i-16} before modulo 2^64 */
        unreduced_W[i].allocate(pb, FMT(this->annotation_prefix, " unreduced_W_%zu", i));

        /* allocate the bit representation of packed_W[i] */
        W_bits[i].allocate(pb, 64, FMT(this->annotation_prefix, " W_bits_%zu", i));

        /* and finally reduce this into packed and bit representations */
        mod_reduce_W[i].reset(new lastbits_gadget<FieldT>(pb, unreduced_W[i], 64+2, packed_W[i], W_bits[i], FMT(this->annotation_prefix, " mod_reduce_W_%zu", i)));
    }
}

template<typename FieldT>
void sha512_message_schedule_gadget<FieldT>::generate_r1cs_constraints()
{
    for (size_t i = 0; i < 16; ++i)
    {
        pack_W[i]->generate_r1cs_constraints(false); // do not enforce bitness here; caller be aware.
    }

    for (size_t i = 16; i < 80; ++i) 
    {
        compute_sigma0[i]->generate_r1cs_constraints();
        compute_sigma1[i]->generate_r1cs_constraints();

        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(1,
                                                             sigma0[i] + sigma1[i] + packed_W[i-16] + packed_W[i-7],
                                                             unreduced_W[i]),
            FMT(this->annotation_prefix, " unreduced_W_%zu", i));

        mod_reduce_W[i]->generate_r1cs_constraints();
    }
}

template<typename FieldT>
void sha512_message_schedule_gadget<FieldT>::generate_r1cs_witness()
{
    for (size_t i = 0; i < 16; ++i)
    {
        pack_W[i]->generate_r1cs_witness_from_bits();
    }

    for (size_t i = 16; i < 80; ++i)
    {
        compute_sigma0[i]->generate_r1cs_witness();
        compute_sigma1[i]->generate_r1cs_witness();

        this->pb.val(unreduced_W[i]) = this->pb.val(sigma0[i]) + this->pb.val(sigma1[i]) + this->pb.val(packed_W[i-16]) + this->pb.val(packed_W[i-7]);
        mod_reduce_W[i]->generate_r1cs_witness();
    }
}

template<typename FieldT>
sha512_round_function_gadget<FieldT>::sha512_round_function_gadget(protoboard<FieldT> &pb,
                                                                   const pb_linear_combination_array<FieldT> &a,
                                                                   const pb_linear_combination_array<FieldT> &b,
                                                                   const pb_linear_combination_array<FieldT> &c,
                                                                   const pb_linear_combination_array<FieldT> &d,
                                                                   const pb_linear_combination_array<FieldT> &e,
                                                                   const pb_linear_combination_array<FieldT> &f,
                                                                   const pb_linear_combination_array<FieldT> &g,
                                                                   const pb_linear_combination_array<FieldT> &h,
                                                                   const pb_variable<FieldT> &W,
                                                                   const long &K,
                                                                   const pb_linear_combination_array<FieldT> &new_a,
                                                                   const pb_linear_combination_array<FieldT> &new_e,
                                                                   const std::string &annotation_prefix) :
    gadget<FieldT>(pb, annotation_prefix),
    a(a),
    b(b),
    c(c),
    d(d),
    e(e),
    f(f),
    g(g),
    h(h),
    W(W),
    K(K),
    new_a(new_a),
    new_e(new_e)
{
    /* compute sigma0 and sigma1 */
    sigma0.allocate(pb, FMT(this->annotation_prefix, " sigma0"));
    sigma1.allocate(pb, FMT(this->annotation_prefix, " sigma1"));
    compute_sigma0.reset(new big_sigma_gadget<FieldT>(pb, a, sigma0, 28, 34, 39, FMT(this->annotation_prefix, " compute_sigma0")));
    compute_sigma1.reset(new big_sigma_gadget<FieldT>(pb, e, sigma1, 14, 18, 41, FMT(this->annotation_prefix, " compute_sigma1")));

    /* compute choice */
    choice.allocate(pb, FMT(this->annotation_prefix, " choice"));
    compute_choice.reset(new choice_gadget<FieldT>(pb, e, f, g, choice, FMT(this->annotation_prefix, " compute_choice")));

    /* compute majority */
    majority.allocate(pb, FMT(this->annotation_prefix, " majority"));
    compute_majority.reset(new majority_gadget<FieldT>(pb, a, b, c, majority, FMT(this->annotation_prefix, " compute_majority")));

    /* pack d */
    packed_d.allocate(pb, FMT(this->annotation_prefix, " packed_d"));
    pack_d.reset(new packing_gadget<FieldT>(pb, d, packed_d, FMT(this->annotation_prefix, " pack_d")));

    /* pack h */
    packed_h.allocate(pb, FMT(this->annotation_prefix, " packed_h"));
    pack_h.reset(new packing_gadget<FieldT>(pb, h, packed_h, FMT(this->annotation_prefix, " pack_h")));

    /* compute the actual results for the round */
    unreduced_new_a.allocate(pb, FMT(this->annotation_prefix, " unreduced_new_a"));
    unreduced_new_e.allocate(pb, FMT(this->annotation_prefix, " unreduced_new_e"));

    packed_new_a.allocate(pb, FMT(this->annotation_prefix, " packed_new_a"));
    packed_new_e.allocate(pb, FMT(this->annotation_prefix, " packed_new_e"));

    mod_reduce_new_a.reset(new lastbits_gadget<FieldT>(pb, unreduced_new_a, 64+3, packed_new_a, new_a, FMT(this->annotation_prefix, " mod_reduce_new_a")));
    mod_reduce_new_e.reset(new lastbits_gadget<FieldT>(pb, unreduced_new_e, 64+3, packed_new_e, new_e, FMT(this->annotation_prefix, " mod_reduce_new_e")));
}

template<typename FieldT>
void sha512_round_function_gadget<FieldT>::generate_r1cs_constraints()
{
    compute_sigma0->generate_r1cs_constraints();
    compute_sigma1->generate_r1cs_constraints();

    compute_choice->generate_r1cs_constraints();
    compute_majority->generate_r1cs_constraints();

    pack_d->generate_r1cs_constraints(false);
    pack_h->generate_r1cs_constraints(false);

    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(1,
                                                         packed_h + sigma1 + choice + K + W + sigma0 + majority,
                                                         unreduced_new_a),
        FMT(this->annotation_prefix, " unreduced_new_a"));

    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(1,
                                                         packed_d + packed_h + sigma1 + choice + K + W,
                                                         unreduced_new_e),
        FMT(this->annotation_prefix, " unreduced_new_e"));

    mod_reduce_new_a->generate_r1cs_constraints();
    mod_reduce_new_e->generate_r1cs_constraints();
}

template<typename FieldT>
void sha512_round_function_gadget<FieldT>::generate_r1cs_witness()
{
    compute_sigma0->generate_r1cs_witness();
    compute_sigma1->generate_r1cs_witness();

    compute_choice->generate_r1cs_witness();
    compute_majority->generate_r1cs_witness();

    pack_d->generate_r1cs_witness_from_bits();
    pack_h->generate_r1cs_witness_from_bits();

    this->pb.val(unreduced_new_a) = this->pb.val(packed_h) + this->pb.val(sigma1) + this->pb.val(choice) + FieldT(K) + this->pb.val(W) + this->pb.val(sigma0) + this->pb.val(majority);
    this->pb.val(unreduced_new_e) = this->pb.val(packed_d) + this->pb.val(packed_h) + this->pb.val(sigma1) + this->pb.val(choice) + FieldT(K) + this->pb.val(W);

    mod_reduce_new_a->generate_r1cs_witness();
    mod_reduce_new_e->generate_r1cs_witness();
}

} // libsnark

#endif // SHA512_COMPONENTS_TCC_
