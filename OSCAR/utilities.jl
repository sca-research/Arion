function legendre_symbol(x, p)
    p = convert(Int, p)
    e = convert(Int, (p - 1) / 2)
    val = x^e
    if typeof(val) == FqFieldElem
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

function generate_field_equations(variables)
    q = size(base_ring(variables[1]))
    field_equations = map(var -> var^q - var, variables)
    return field_equations
end
