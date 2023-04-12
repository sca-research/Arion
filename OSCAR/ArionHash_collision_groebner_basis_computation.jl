using Oscar
include("ArionHash_polynomial_model.jl")
include("zero_dimensional_ideals.jl")
include("utilities.jl")

parameters = [
# [field, branches, rounds, capacity, d_2]
[1013, 2, 1, 1, 7],
[1013, 3, 1, 2, 7],
[1013, 4, 1, 3, 7],
[1013, 5, 1, 4, 7],

[10007, 2, 1, 1, 7],
[10007, 3, 1, 2, 7],
[10007, 4, 1, 3, 7],
[10007, 5, 1, 4, 7],

[1033, 2, 1, 1, 7],
[1033, 3, 1, 2, 7],
[1033, 4, 1, 3, 7],
[1033, 5, 1, 4, 7],

[15013, 2, 1, 1, 7],
[15013, 3, 1, 2, 7],
[15013, 4, 1, 3, 7],
[15013, 5, 1, 4, 7],

[1013, 2, 1, 1, 257],
[1013, 3, 1, 2, 257],

[10007, 2, 1, 1, 257],
[10007, 3, 1, 2, 257],

[1033, 2, 1, 1, 257],
[1033, 3, 1, 2, 257],

[15013, 2, 1, 1, 257],
[15013, 3, 1, 2, 257],
]

for param in parameters
    arion_hash = ArionHash_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], capacity=param[4], d_2=param[5])
    polys = generate_ArionHash_collision_polynomials(arion_hash=arion_hash)
    variables = gens(parent(polys[1]))
    lms = get_leading_monomials(polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", get_maximum_degree(lms))
    I = ideal(polys)
    for i in 0:param[2] - 2
        println("Guess for variable: ", variables[length(variables) - arion_hash.rounds - i])
        el = rand(arion_hash.field)
        println("Random guess: ", el)
        J = ideal([variables[length(variables) - arion_hash.rounds] - el])
        gb = f4(I + J, initial_hts=17, nr_thrds=Threads.nthreads(), max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
        println("Size of Gröbner basis: ", length(gb))
        basis = vector_space_basis(gb)
        println("Size of vector space basis: ", length(basis))
        println("Leading monomials of Gröbner basis: ", get_leading_monomials(gb))
        println("")
    end
end
