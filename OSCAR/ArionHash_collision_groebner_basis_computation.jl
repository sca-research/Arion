using Oscar
include("ArionHash_polynomial_model.jl")
include("zero_dimensional_ideals.jl")
include("utilities.jl")

parameters = [
# [field, branches, rounds, capacity, d_2]
[1013, 2, 1, 1, 3],
[1013, 2, 2, 1, 3],

[1013, 3, 1, 2, 3],
[1013, 4, 1, 3, 3],

[10007, 2, 1, 1, 3],
[10007, 2, 2, 1, 3],

[10007, 3, 1, 2, 3],
[10007, 4, 1, 3, 3],

[1013, 2, 1, 1, 7],
[1013, 3, 1, 2, 7],
[1013, 4, 1, 3, 7],

[10007, 2, 1, 1, 7],
[10007, 3, 1, 2, 7],
[10007, 4, 1, 3, 7],

[1033, 2, 1, 1, 5],
[1033, 3, 1, 2, 5],

[15013, 2, 1, 1, 5],
[15013, 3, 1, 2, 5],

[1033, 2, 1, 1, 7],
[1033, 3, 1, 2, 7],

[15013, 2, 1, 1, 7],
[15013, 3, 1, 2, 7],
]

println("Matrix Type 1")
println("\n")
for param in parameters
    arion_hash = ArionHash_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], capacity=param[4], d_2=param[5], matrix_type=1)
    polys = generate_ArionHash_collision_polynomials(arion_hash=arion_hash)
    variables = gens(parent(polys[1]))
    lms = map(poly -> leading_monomial(poly), polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", map(mon -> total_degree(mon), lms))
    I = ideal(polys)
    for i in 0:param[2] - 2
        println("Guess for variable: ", variables[length(variables) - arion_hash.rounds - i])
        el = rand(arion_hash.field)
        println("Random guess: ", el)
        J = ideal([variables[length(variables) - arion_hash.rounds] - el])
        gb = groebner_basis_f4(I + J, initial_hts=17, nr_thrds=Threads.nthreads(), max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
        println("Size of Groebner basis: ", length(gb))
        basis = vector_space_basis(gb)
        println("Size of vector space basis: ", length(basis))
        println("Leading monomials of Gröbner basis: ", map(poly -> leading_monomial(poly), gb))
        println("")
    end
end

println("Matrix Type 2")
println("\n")
for param in parameters
    arion_hash = ArionHash_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], capacity=param[4], d_2=param[5], matrix_type=2)
    polys = generate_ArionHash_collision_polynomials(arion_hash=arion_hash)
    variables = gens(parent(polys[1]))
    lms = map(poly -> leading_monomial(poly), polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", map(mon -> total_degree(mon), lms))
    I = ideal(polys)
    for i in 0:param[2] - 2
        println("Guess for variable: ", variables[length(variables) - arion_hash.rounds - i])
        el = rand(arion_hash.field)
        println("Random guess: ", el)
        J = ideal([variables[length(variables) - arion_hash.rounds] - el])
        gb = groebner_basis_f4(I + J, initial_hts=17, nr_thrds=Threads.nthreads(), max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
        println("Size of Groebner basis: ", length(gb))
        basis = vector_space_basis(gb)
        println("Size of vector space basis: ", length(basis))
        println("Leading monomials of Gröbner basis: ", map(poly -> leading_monomial(poly), gb))
        println("")
    end
end
