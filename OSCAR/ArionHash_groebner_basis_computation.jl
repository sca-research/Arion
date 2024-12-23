using Oscar
include("ArionHash_polynomial_model.jl")
include("zero_dimensional_ideals.jl")
include("utilities.jl")

nr_thrds = 16

parameters = [
# [field, branches, rounds, capacity, d_2]
[1013, 2, 1, 1, 7],
[1013, 2, 2, 1, 7],
[1013, 2, 3, 1, 7],

[1013, 3, 1, 2, 7],
[1013, 3, 2, 2, 7],

[10007, 2, 1, 1, 7],
[10007, 2, 2, 1, 7],
[10007, 2, 3, 1, 7],

[10007, 3, 1, 2, 7],
[10007, 3, 2, 2, 7],

[1033, 2, 1, 1, 7],
[1033, 2, 2, 1, 7],
[1033, 2, 3, 1, 7],

[1033, 3, 1, 2, 7],
[1033, 3, 2, 2, 7],

[15013, 2, 1, 1, 7],
[15013, 2, 2, 1, 7],
[15013, 2, 3, 1, 7],

[15013, 3, 1, 2, 7],
[15013, 3, 2, 2, 7],
]

println("Matrix Type 1")
println("\n")
for param in parameters
    arion_hash = ArionHash_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], capacity=param[4], d_2=param[5], matrix_type=1)
    polys = generate_ArionHash_polynomials(arion_hash=arion_hash)
    lms = map(poly -> leading_monomial(poly), polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", map(mon -> total_degree(mon), lms))
    I = ideal(polys)
    gb = groebner_basis_f4(I, initial_hts=17, nr_thrds=nr_thrds, max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
    println("Size of Gröbner basis: ", length(gb))
    basis = vector_space_basis(gb)
    println("Size of vector space basis: ", length(basis))
    println("Leading monomials of Gröbner basis: ", map(poly -> leading_monomial(poly), gb))
    println("")
end

println("Matrix Type 2")
println("\n")
for param in parameters
    arion_hash = ArionHash_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], capacity=param[4], d_2=param[5], matrix_type=2)
    polys = generate_ArionHash_polynomials(arion_hash=arion_hash)
    lms = map(poly -> leading_monomial(poly), polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", map(mon -> total_degree(mon), lms))
    I = ideal(polys)
    gb = groebner_basis_f4(I, initial_hts=17, nr_thrds=nr_thrds, max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
    println("Size of Gröbner basis: ", length(gb))
    basis = vector_space_basis(gb)
    println("Size of vector space basis: ", length(basis))
    println("Leading monomials of Gröbner basis: ", map(poly -> leading_monomial(poly), gb))
    println("")
end
