function legendre_symbol(x, p)
    p = convert(Int, p)
    e = convert(Int, (p - 1) / 2)
    val = x^e
    if typeof(val) == gfp_elem
        if val == 0
            return 0
        elseif val == 1
            return 1
        else
            return -1
        end
    end
    val = val % p
    if val == 0.0
        return 0
    elseif val == 1.0
        return 1
    else
        return -1
    end
end

function get_leading_monomials(polys)
    lms = Vector{typeof(polys[1])}()
    try
        lms = copy(polys)
    catch
        lms = copy(gens(polys))
    end
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

function append_polynomial_matrix_to_vector(poly_vector, poly_matrix)
    for i in 1:nrows(poly_matrix)
        temp_poly = poly_matrix[i, 1]
        push!(poly_vector, temp_poly)
    end
    return poly_vector
end

function generate_field_equations(variables)
    q = size(base_ring(variables[1]))
    field_equations = Vector{typeof(variables[1])}()
    for i in 1:length(variables)
        push!(field_equations, variables[i]^q - variables[i])
    end
    return field_equations
end

function int_vector_to_field_matrix(v_in, field)
    v_out = zero_matrix(field, length(v_in), 1)
    for i in 1:length(v_in)
        try
            v_out[i, 1] = field(v_in[i])
        catch
            v_out[i, 1] = field(v_in[i, 1])
        end
    end
    return v_out
end
