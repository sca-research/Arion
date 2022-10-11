using Oscar
include("Arion.jl")

function append_polynomial_matrix_to_vector(poly_vector, poly_matrix)
    for i in 1:nrows(poly_matrix)
        temp_poly = poly_matrix[i, 1]
        push!(poly_vector, temp_poly)
    end
    return poly_vector
end

function generate_variables(branches, rounds)
    variables = String[]
    push!(variables, "z_" * string(1))
    for i in 1:(rounds - 1)
        for j in 1:branches
            push!(variables, "x_" * string(j) * "__" * string(i))
        end
        push!(variables, "z_" * string(i + 1))
    end
    for i in 1:branches
        push!(variables, "y_" * string(i))
    end
    return variables
end

function generate_field_equations(variables)
    q = size(base_ring(variables[1]))
    field_equations = Vector{typeof(variables[1])}()
    for i in 1:length(variables)
        push!(field_equations, variables[i]^q - variables[i])
    end
    return field_equations
end

function GTDS_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
    v_out = zero_matrix(base_ring(v_in), branches, 1)
    v_out[branches, 1] = v_in[branches, 1]
    sigma = v_n_std + v_out[branches, 1]
    for i in branches - 1:-1:1
        v_out[i, 1] = v_in[i, 1]^d
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i,:], constants_h[i,:])
        v_out[i, 1] *= g_i
        v_out[i, 1] += h_i
        sigma += v_in[i, 1] + v_out[i, 1]
    end
    return v_out
end

function round_function_polynomial_model(v_in, v_n_std, key, branches, d, constants_g, constants_h, constants_aff, mat)
    v_out = GTDS_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
    v_out = affine_layer(v_out, constants_aff, mat)
    v_out = key_addition(v_out, key)
    return v_out
end

function generate_Arion_polynomials(;arion=Arion_constructor(),
                                     plain=nothing,
                                     cipher=nothing,
                                     termorder="degrevlex",
                                     field_equations=false,
                                     naive_model=false)
    print_key = false
    if plain == nothing
        plain = zero_matrix(arion.field, arion.branches, 1)
        for i in 1:arion.branches
            plain[i, 1] = rand(arion.field)
        end
        if cipher == nothing
            print_key = true
            key = zero_matrix(arion.field, arion.branches, 1)
            for i in 1:arion.branches
                key[i, 1] = rand(arion.field)
            end
            cipher = encrypt(plain, key, arion)
        end
    end

    println("Plain text: ", plain)
    if print_key
        println("Key: ", key)
    end
    println("Cipher text: ", cipher)
    println("Term order: ", termorder)

    if termorder == "degrevlex"
        P, variables = PolynomialRing(arion.field, generate_variables(arion.branches, arion.rounds), ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = PolynomialRing(arion.field, generate_variables(arion.branches, arion.rounds), ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end
    
    polynomials = Vector{typeof(variables[1])}()
    key_variables = matrix(variables[arion.branches * (arion.rounds - 1) + arion.rounds + 1:arion.branches * arion.rounds + arion.rounds])
    if arion.rounds > 1
        index = 1
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[index]
        next_state = matrix(variables[index + 1:arion.branches + 1])
        polynomials = append_polynomial_matrix_to_vector(polynomials, round_function_polynomial_model(current_state,
                                                                                                      tmp,
                                                                                                      key_variables,
                                                                                                      arion.branches,
                                                                                                      arion.d_1,
                                                                                                      arion.constants_g[1:arion.branches,:],
                                                                                                      arion.constants_h[1:arion.branches,:],
                                                                                                      arion.constants_aff[1,:],
                                                                                                      arion.matrix) - next_state)
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion.d_2)
        end
        index += 1
        for r in 1:arion.rounds - 2
            current_state = matrix(variables[(r - 1) * arion.branches + index:r * arion.branches + index - 1])
            tmp = current_state[arion.branches, 1]
            current_state[arion.branches, 1] = variables[r * arion.branches + index]
            next_state = matrix(variables[r * arion.branches + index + 1:(r + 1) * arion.branches + index])
            polynomials = append_polynomial_matrix_to_vector(polynomials,
                                                             round_function_polynomial_model(current_state,
                                                                                             tmp,
                                                                                             key_variables,
                                                                                             arion.branches,
                                                                                             arion.d_1,
                                                                                             arion.constants_g[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                                             arion.constants_h[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                                             arion.constants_aff[r + 1,:],
                                                                                             arion.matrix) - next_state)
            if naive_model
                push!(polynomials, tmp^arion.e_2 - variables[r * arion.branches + index])
            else
                push!(polynomials, tmp - variables[r * arion.branches + index]^arion.d_2)
            end
            index += 1
        end
        current_state = matrix(variables[(arion.rounds - 2) * arion.branches + index:(arion.rounds - 1) * arion.branches + index - 1])
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[(arion.rounds - 1) * arion.branches + index]
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_polynomial_model(current_state,
                                                                                         tmp,
                                                                                         key_variables,
                                                                                         arion.branches,
                                                                                         arion.d_1,
                                                                                         arion.constants_g[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                                                         arion.constants_h[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                                                         arion.constants_aff[arion.rounds,:],
                                                                                         arion.matrix) - cipher)
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[(arion.rounds - 1) * arion.branches + index])
        else
            push!(polynomials, tmp - variables[(arion.rounds - 1) * arion.branches + index]^arion.d_2)
        end
    else
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[1]
        polynomials = append_polynomial_matrix_to_vector(polynomials, round_function_polynomial_model(current_state,
                                                                                                      tmp,
                                                                                                      key_variables,
                                                                                                      arion.branches,
                                                                                                      arion.d_1,
                                                                                                      arion.constants_g,
                                                                                                      arion.constants_h,
                                                                                                      arion.constants_aff[1,:],
                                                                                                      arion.matrix) - cipher)
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[1])
        else
            push!(polynomials, tmp - variables[1]^arion.d_2)
        end
    end
    
    if field_equations
        fes = generate_field_equations(variables)
        polynomials = append_polynomial_matrix_to_vector(polynomials, matrix(fes))
    end
    return polynomials
end
