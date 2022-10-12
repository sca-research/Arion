using Oscar
include("Arion_polynomial_model.jl")
include("zero_dimensional_ideals.jl")

function get_leading_monomials(polys)
    lms = copy(polys)
    for i in 1:length(polys)
        lms[i] = leading_monomial(lms[i])
    end
    return lms
end

function get_maximum_degree(polys)
    max_deg = 0
    for i in 1:length(polys)
        d = total_degree(polys[i])
        if d > max_deg
            max_deg = d
        end
    end
    return max_deg
end

parameters = [
# [field, branches, rounds, d_2]
[1013, 2, 1, 7],
[1013, 2, 1, 257],
[1013, 3, 1, 7],
[1013, 4, 1, 7],

[10007, 2, 1, 7],
[10007, 2, 1, 257],
[10007, 3, 1, 7],
[10007, 4, 1, 7],

[1033, 2, 1, 7],
[1033, 2, 1, 257],
[1033, 3, 1, 7],

[15013, 2, 1, 7],
[15013, 2, 1, 257],
[15013, 3, 1, 7],
]

for param in parameters
    arion = Arion_constructor(field=GF(param[1]), branches=param[2], rounds=param[3], d_2=param[4])
    polys = generate_Arion_polynomials(arion=arion)
    lms = get_leading_monomials(polys)
    println("Leading monomials: ", lms)
    println("Maximum degree: ", get_maximum_degree(lms))
    I = ideal(polys)
    gb = f4(I, initial_hts=17, nr_thrds=Threads.nthreads(), max_nr_pairs=0, la_option=2, eliminate=0, complete_reduction=true, info_level=2)
    println("Size of Gröbner basis: ", length(gb))
    basis = vector_space_basis(gb)
    println("Size of vector space basis: ", length(basis))
    println("Leading monomials of Gröbner basis: ", get_leading_monomials(gb))
    println("")
end
