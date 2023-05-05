using Oscar
include("zero_dimensional_ideals.jl")

function change_of_term_order_for_Arion_generic(P, branches, rounds)
    variables = gens(P)
    variables_P_1 = variables[1:branches * (rounds - 1)]
    variables_P_2 = variables[branches * (rounds - 1) + 1:branches * rounds]
    variables_P_3 = variables[branches * rounds + 1:length(variables)]

    variables_Q = String[]
    for var in [variables_P_2; variables_P_1; variables_P_3]
        push!(variables_Q, string(var))
    end
    Q, variables_Q = PolynomialRing(base_ring(P), variables_Q, ordering=ordering(P))
    variables = gens(Q)
    variables_Q_1 = variables[1:branches]
    variables_Q_2 = variables[branches + 1:branches * rounds]
    variables_Q_3 = variables[branches * rounds + 1:length(variables)]

    F = hom(P, Q, [variables_Q_2; variables_Q_1; variables_Q_3])
    G = hom(Q, P, [variables_P_2; variables_P_1; variables_P_3])

    return F, G
end

function compute_Arion_decomposition_ideal_quotient_space_dimension(decomposition_leading_terms, decomposition_variables)
    I_leading_terms = copy(decomposition_leading_terms)
    poly_ring = parent(decomposition_leading_terms[1])
    for var in gens(poly_ring)
        if !(var in decomposition_variables)
            push!(I_leading_terms, var)
        end
    end
    I_leading_terms = ideal(poly_ring, I_leading_terms)
    if dim(I_leading_terms) == 0
        return length(vector_space_basis(gens(I_leading_terms), sort_basis=false))
    else
        println("Ideal is not zero-dimensional!")
        return -1
    end
end

function Arion_groebner_basis_compution(arion, arion_polynomial_system; nr_thrds=1, info_level=0)
    P = parent(arion_polynomial_system[1])
    variables = gens(P)
    F, G = change_of_term_order_for_Arion_generic(P, arion.branches, arion.rounds);
    Q = parent(F(variables[1]));
    variables_Q = gens(Q);

    variable_parition =  Vector{typeof([variables[1]])}()
    decomposition_ideals = Vector{typeof(ideal(variables))}()
    decomposition_groebner_bases = Vector{typeof(groebner_basis_f4(ideal(variables)))}()
    decomposition_quotient_space_dimensions = Vector{Int128}()
    linear_polys =  Vector{typeof(variables[1])}()

    r = 0
    println("Round: ", r)
    X_part = Vector{typeof(variables[1])}()
    I_decomp = Vector{typeof(variables[1])}()

    for r in 1:arion.rounds
        push!(X_part, variables_Q[arion.rounds * arion.branches + r])
        push!(I_decomp, arion_polynomial_system[r * (arion.branches + 1)])
    end
    push!(variable_parition, X_part)
    println("Variable partition:")
    for var in X_part
        println(var)
    end

    I_decomp = F(ideal(P, I_decomp))
    push!(decomposition_ideals, I_decomp)
    gb_decomp = groebner_basis_f4(I_decomp)
    push!(decomposition_groebner_bases, gb_decomp)
    dim_decomp = BigInt(arion.d_2)^arion.rounds
    push!(decomposition_quotient_space_dimensions, dim_decomp)
    println("Leading monomials:")
    for mon in map(poly -> leading_monomial(poly), gb_decomp)
        println(mon)
    end
    println("Quotient space dimension: ", dim_decomp)
    println("\n")

    r = 1
    println("Round: ", r)
    X_part = variables_Q[(r - 1) * arion.branches + 1:r * arion.branches]
    push!(variable_parition, X_part)
    println("Variable partition:")
    for var in X_part
        println(var)
    end

    I_decomp = F(ideal(P, arion_polynomial_system[1:r * arion.branches]))
    push!(decomposition_ideals, I_decomp)

    gb_decomp = groebner_basis_f4(I_decomp, nr_thrds=nr_thrds, info_level=info_level)
    push!(decomposition_groebner_bases, gb_decomp)
    println("Leading monomials:")
    gb_decomp_leading_monomials = map(poly -> leading_monomial(poly), gb_decomp)
    for mon in gb_decomp_leading_monomials
        println(mon)
    end
    dim_decomp = compute_Arion_decomposition_ideal_quotient_space_dimension(gb_decomp_leading_monomials, X_part)
    push!(decomposition_quotient_space_dimensions, dim_decomp)
    println("Quotient space dimension: ", dim_decomp)
    println("\n")

    linear_polys = [linear_polys; collect(Iterators.filter(poly -> total_degree(poly) == 1, gb_decomp))]
    I_forward = ideal(Q, linear_polys)

    for r in 2:arion.rounds
        println("Round: ", r)
        X_part = variables_Q[(r - 1) * arion.branches + 1:r * arion.branches]
        push!(variable_parition, X_part)
        println("Variable partition:")
        for var in X_part
            println(var)
        end

        I_decomp = F(ideal(P, arion_polynomial_system[(r - 1) * (arion.branches + 1) + 1:(r - 1) * (arion.branches + 1) + arion.branches])) + I_forward
        push!(decomposition_ideals, I_decomp)

        gb_decomp = groebner_basis_f4(I_decomp, nr_thrds=nr_thrds, info_level=info_level)
        push!(decomposition_groebner_bases, gb_decomp)
        println("Leading monomials:")
        gb_decomp_leading_monomials = map(poly -> leading_monomial(poly), gb_decomp)
        for mon in gb_decomp_leading_monomials
            println(mon)
        end
        dim_decomp = compute_Arion_decomposition_ideal_quotient_space_dimension(gb_decomp_leading_monomials, X_part)
        push!(decomposition_quotient_space_dimensions, dim_decomp)
        println("Quotient space dimension: ", dim_decomp)
        println("\n")

        linear_polys = [linear_polys; collect(Iterators.filter(poly -> total_degree(poly) == 1, gb_decomp))]
        I_forward = ideal(Q, linear_polys)
    end

    gb_Arion = Vector{typeof(variables[1])}()
    for (gb, X_part) in zip(decomposition_groebner_bases, variable_parition)
        for poly in gens(gb)
            for var in X_part
                if divides(leading_monomial(poly), var)[1]
                    push!(gb_Arion, poly)
                    break
                end
            end
        end
    end

    if !(-1 in decomposition_quotient_space_dimensions)
        D::BigInt = prod(decomposition_quotient_space_dimensions)
    else
        D = -1
    end

    return gb_Arion, D
end
