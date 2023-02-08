using Oscar
using Folds

function push_vec(vec_1, vec_2)
    for i in 1:length(vec_2)
        push!(vec_1, vec_2[i])
    end
    return vec_1
end

function divisible_by_leading_monomial(monom, leading_monoms)
    poly_ring = parent(monom)
    division = unique(Folds.map(x -> div(monom, x), leading_monoms))
    if division == [poly_ring(0)]
        return monom
    else
        return poly_ring(1)
    end
end

function vector_space_basis(gb; sort_basis=true)
    poly_ring = parent(gb[1])
    variables = gens(poly_ring)
    lms = Vector{typeof(gb[1])}()
    try
        lms = copy(gb)
    catch
        lms = copy(gens(gb))
    end
    for i in 1:length(gb)
        lms[i] = leading_monomial(lms[i])
    end
    basis = [poly_ring(1)]
    basis_found = false
    prevous_new_basis_elements = copy(basis)
    while !basis_found
        cart_iter = Iterators.product(variables, prevous_new_basis_elements)
        cart_prod = vec(Folds.map(x -> x[1] * x[2], cart_iter))
        new_basis_elements = unique(Folds.map(x -> divisible_by_leading_monomial(x, lms), cart_prod))
        deleteat!(new_basis_elements, findall(x -> x in basis, new_basis_elements))
        if new_basis_elements == []
            basis_found = true
        else
            basis = push_vec(basis, new_basis_elements)
            prevous_new_basis_elements = copy(new_basis_elements)
            new_basis_elements = Vector{typeof(variables[1])}()
        end
    end
    basis = unique(basis)
    if sort_basis
        return sort(basis, rev=true)
    else
        return basis
    end
end
