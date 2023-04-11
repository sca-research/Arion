using Oscar
include("ArionHash.jl")
include("utilities.jl")

function generate_ArionHash_variables(branches, rounds, rate)
    variables = String[]
    for i in 1:rate
        push!(variables, "x_in__" * string(i))
    end
    push!(variables, "z_" * string(1))
    for i in 1:rounds - 1
        for j in 1:branches
            push!(variables, "x_" * string(j) * "__" * string(i))
        end
        push!(variables, "z_" * string(i + 1))
    end
    for i in 2:branches
        push!(variables, "x_out__" * string(i))
    end
    return variables
end

function GTDS_ArionHash_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
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

function round_function_ArionHash_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h, constants_aff, mat)
    v_out = GTDS_ArionHash_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
    v_out = affine_layer(v_out, constants_aff, mat)
    return v_out
end

function generate_ArionHash_polynomials(;arion_hash=ArionHash_constructor(),
                                         hash_val=nothing,
                                         termorder="degrevlex",
                                         field_equations=false,
                                         naive_model=false)
    print_plain = false
    if hash_val == nothing
        print_plain = true
        plain = zero_matrix(arion_hash.field, arion_hash.rate, 1)
        for i in 1:arion_hash.rate
            plain[i, 1] = rand(arion_hash.field)
        end
        hash_val = hash(plain, arion_hash)
    else
        hash_val = arion_hash.field(hash_val)
    end

    if print_plain
        println("Plain text: ", plain)
    end
    println("Hash value: ", hash_val)

    variables = generate_ArionHash_variables(arion_hash.branches, arion_hash.rounds, arion_hash.rate)
    if termorder == "degrevlex"
        P, variables = PolynomialRing(arion_hash.field, variables, ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = PolynomialRing(arion_hash.field, variables, ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end
    
    polynomials = Vector{typeof(variables[1])}()
    current_state = zero_matrix(P, arion_hash.branches, 1)

    for i in 1:arion_hash.branches
        if i <= arion_hash.rate
            current_state[i, 1] += variables[i]
        else
            current_state[i, 1] += arion_hash.initial_value[i, 1]
        end
    end
    index = arion_hash.rate + 1
    current_state = arion_hash.matrix * current_state
    tmp = current_state[arion_hash.branches, 1]
    current_state[arion_hash.branches, 1] = variables[index]
    if arion_hash.rounds > 1
        next_state = matrix(variables[index + 1:index + arion_hash.branches])
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_ArionHash_polynomial_model(current_state,
                                                                                                   tmp,
                                                                                                   arion_hash.branches,
                                                                                                   arion_hash.d_1,
                                                                                                   arion_hash.constants_g[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_h[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_aff[1,:],
                                                                                                   arion_hash.matrix) - next_state)
        if naive_model
            push!(polynomials, tmp^arion_hash.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion_hash.d_2)
        end
        for r in 1:arion_hash.rounds - 2
            current_state = matrix(variables[(r - 1) * arion_hash.branches + index + 1:r * arion_hash.branches + index])
            tmp = current_state[arion_hash.branches, 1]
            index += 1
            current_state[arion_hash.branches, 1] = variables[r * arion_hash.branches + index]
            next_state = matrix(variables[r * arion_hash.branches + index + 1:(r + 1) * arion_hash.branches + index])
            polynomials = append_polynomial_matrix_to_vector(polynomials,
                                                             round_function_ArionHash_polynomial_model(current_state,
                                                                                                       tmp,
                                                                                                       arion_hash.branches,
                                                                                                       arion_hash.d_1,
                                                                                                       arion_hash.constants_g[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                       arion_hash.constants_h[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                       arion_hash.constants_aff[r + 1,:],
                                                                                                       arion_hash.matrix) - next_state)
            if naive_model
                push!(polynomials, tmp^arion_hash.e_2 - variables[r * arion_hash.branches + index])
            else
                push!(polynomials, tmp - variables[r * arion_hash.branches + index]^arion_hash.d_2)
            end
        end
        current_state = matrix(variables[(arion_hash.rounds - 2) * arion_hash.branches + index + 1:(arion_hash.rounds - 1) * arion_hash.branches + index])
        tmp = current_state[arion_hash.branches, 1]
        index += 1
        current_state[arion_hash.branches, 1] = variables[(arion_hash.rounds - 1) * arion_hash.branches + index]
        next_state = zero_matrix(P, arion_hash.branches, 1)
        next_state[1, 1] = hash_val
        for i in 2:arion_hash.branches
            next_state[i, 1] = variables[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
        end
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_ArionHash_polynomial_model(current_state,
                                                                                                   tmp,
                                                                                                   arion_hash.branches,
                                                                                                   arion_hash.d_1,
                                                                                                   arion_hash.constants_g[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                                                                   arion_hash.constants_h[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                                                                   arion_hash.constants_aff[arion_hash.rounds,:],
                                                                                                   arion_hash.matrix) - next_state)
        if naive_model
            push!(polynomials, tmp^arion_hash.e_2 - variables[(arion_hash.rounds - 1) * arion_hash.branches + index])
        else
            push!(polynomials, tmp - variables[(arion_hash.rounds - 1) * arion_hash.branches + index]^arion_hash.d_2)
        end
    else
        next_state = zero_matrix(P, arion_hash.branches, 1)
        next_state[1, 1] = hash_val
        for i in 2:arion_hash.branches
            next_state[i, 1] = variables[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
        end
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_ArionHash_polynomial_model(current_state,
                                                                                                   tmp,
                                                                                                   arion_hash.branches,
                                                                                                   arion_hash.d_1,
                                                                                                   arion_hash.constants_g,
                                                                                                   arion_hash.constants_h,
                                                                                                   arion_hash.constants_aff[1,:],
                                                                                                   arion_hash.matrix) - next_state)
        if naive_model
            push!(polynomials, tmp^arion_hash.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion_hash.d_2)
        end
    end

    if field_equations
        fes = generate_field_equations(variables)
        polynomials = append_polynomial_matrix_to_vector(polynomials, matrix(fes))
    end
    return polynomials
end

function generate_ArionHash_collision_variables(branches, rounds, rate)
    variables_1 = String[]
    for i in 1:rate
        push!(variables_1, "x_in_1__" * string(i))
    end
    push!(variables_1, "z_1_" * string(1))
    for i in 1:rounds - 1
        for j in 1:branches
            push!(variables_1, "x_1_" * string(j) * "__" * string(i))
        end
        push!(variables_1, "z_1_" * string(i + 1))
    end
    for i in 2:branches
        push!(variables_1, "x_out_1__" * string(i))
    end
    variables_2 = String[]
    for i in 1:rate
        push!(variables_2, "x_in_2__" * string(i))
    end
    push!(variables_2, "z_2_" * string(1))
    for i in 1:rounds - 1
        for j in 1:branches
            push!(variables_2, "x_2_" * string(j) * "__" * string(i))
        end
        push!(variables_2, "z_2_" * string(i + 1))
    end
    for i in 2:branches
        push!(variables_2, "x_out_2__" * string(i))
    end
    return [variables_1; variables_2]
end

function generate_ArionHash_collision_polynomials(;arion_hash=ArionHash_constructor(),
                                                    termorder="degrevlex",
                                                    field_equations=false,
                                                    naive_model=false)

    variables = generate_ArionHash_collision_variables(arion_hash.branches, arion_hash.rounds, arion_hash.rate)
    if termorder == "degrevlex"
        P, variables = PolynomialRing(arion_hash.field, variables, ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = PolynomialRing(arion_hash.field, variables, ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end

    N = arion_hash.rate + (arion_hash.branches + 1) * (arion_hash.rounds - 1) + arion_hash.branches
    variables_1 = variables[1:N]
    variables_2 = variables[N + 1:2 * N]

    polynomials = Vector{typeof(variables[1])}()
    current_state_1 = zero_matrix(P, arion_hash.branches, 1)
    current_state_2 = zero_matrix(P, arion_hash.branches, 1)

    for i in 1:arion_hash.branches
        if i <= arion_hash.rate
            current_state_1[i, 1] += variables_1[i]
            current_state_2[i, 1] += variables_2[i]
        else
            current_state_1[i, 1] += arion_hash.initial_value[i, 1]
            current_state_2[i, 1] += arion_hash.initial_value[i, 1]
        end
    end

    index = arion_hash.rate + 1
    current_state_1 = arion_hash.matrix * current_state_1
    current_state_2 = arion_hash.matrix * current_state_2
    tmp_1 = current_state_2[arion_hash.branches, 1]
    tmp_2 = current_state_2[arion_hash.branches, 1]
    current_state_1[arion_hash.branches, 1] = variables_1[index]
    current_state_2[arion_hash.branches, 1] = variables_2[index]
    if arion_hash.rounds > 1
        next_state_1 = matrix(variables_1[index + 1:index + arion_hash.branches])
        next_state_2 = matrix(variables_2[index + 1:index + arion_hash.branches])
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_ArionHash_polynomial_model(current_state_1,
                                                                                                   tmp_1,
                                                                                                   arion_hash.branches,
                                                                                                   arion_hash.d_1,
                                                                                                   arion_hash.constants_g[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_h[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_aff[1,:],
                                                                                                   arion_hash.matrix) - next_state_1)
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_ArionHash_polynomial_model(current_state_2,
                                                                                                   tmp_2,
                                                                                                   arion_hash.branches,
                                                                                                   arion_hash.d_1,
                                                                                                   arion_hash.constants_g[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_h[1:arion_hash.branches,:],
                                                                                                   arion_hash.constants_aff[1,:],
                                                                                                   arion_hash.matrix) - next_state_2)
        if naive_model
            push!(polynomials, tmp_1^arion_hash.e_2 - variables_1[index])
            push!(polynomials, tmp_2^arion_hash.e_2 - variables_2[index])
        else
            push!(polynomials, tmp_1 - variables_1[index]^arion_hash.d_2)
            push!(polynomials, tmp_2 - variables_2[index]^arion_hash.d_2)
        end
        for r in 1:arion_hash.rounds - 2
            current_state_1 = matrix(variables_1[(r - 1) * arion_hash.branches + index + 1:r * arion_hash.branches + index])
            current_state_2 = matrix(variables_2[(r - 1) * arion_hash.branches + index + 1:r * arion_hash.branches + index])
            tmp_1 = current_state_1[arion_hash.branches, 1]
            tmp_2 = current_state_2[arion_hash.branches, 1]
            index += 1
            current_state_1[arion_hash.branches, 1] = variables_1[r * arion_hash.branches + index]
            current_state_2[arion_hash.branches, 1] = variables_2[r * arion_hash.branches + index]
            next_state_1 = matrix(variables_1[r * arion_hash.branches + index + 1:(r + 1) * arion_hash.branches + index])
            next_state_2 = matrix(variables_2[r * arion_hash.branches + index + 1:(r + 1) * arion_hash.branches + index])
            polynomials = append_polynomial_matrix_to_vector(polynomials,
                                                             round_function_ArionHash_polynomial_model(current_state_1,
                                                                                                      tmp_1,
                                                                                                      arion_hash.branches,
                                                                                                      arion_hash.d_1,
                                                                                                      arion_hash.constants_g[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                      arion_hash.constants_h[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                      arion_hash.constants_aff[r + 1,:],
                                                                                                      arion_hash.matrix) - next_state_1)
            polynomials = append_polynomial_matrix_to_vector(polynomials,
                                                             round_function_ArionHash_polynomial_model(current_state_2,
                                                                                                      tmp_2,
                                                                                                      arion_hash.branches,
                                                                                                      arion_hash.d_1,
                                                                                                      arion_hash.constants_g[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                      arion_hash.constants_h[(arion_hash.branches - 1) * r + 1:(arion_hash.branches - 1) * (r + 1),:],
                                                                                                      arion_hash.constants_aff[r + 1,:],
                                                                                                      arion_hash.matrix) - next_state_2)
            if naive_model
                push!(polynomials, tmp_1^arion_hash.e_2 - variables_1[index])
                push!(polynomials, tmp_2^arion_hash.e_2 - variables_2[index])
            else
                push!(polynomials, tmp_1 - variables_1[index]^arion_hash.d_2)
                push!(polynomials, tmp_2 - variables_2[index]^arion_hash.d_2)
            end
        end
        current_state_1 = matrix(variables_1[(arion_hash.rounds - 2) * arion_hash.branches + index + 1:(arion_hash.rounds - 1) * arion_hash.branches + index])
        current_state_2 = matrix(variables_2[(arion_hash.rounds - 2) * arion_hash.branches + index + 1:(arion_hash.rounds - 1) * arion_hash.branches + index])
        tmp_1 = current_state_1[arion_hash.branches, 1]
        tmp_2 = current_state_2[arion_hash.branches, 1]
        index += 1
        current_state_1[arion_hash.branches, 1] = variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index]
        current_state_2[arion_hash.branches, 1] = variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index]
        next_state_1 = zero_matrix(P, arion_hash.branches, 1)
        next_state_2 = zero_matrix(P, arion_hash.branches, 1)
        next_state_1[1, 1] = zero(arion_hash.field)
        next_state_2[1, 1] = zero(arion_hash.field)
        for i in 2:arion_hash.branches
            next_state_1[i, 1] = variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
            next_state_2[i, 1] = variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
        end
        polys_1 = round_function_ArionHash_polynomial_model(current_state_1,
                                                            tmp_1,
                                                            arion_hash.branches,
                                                            arion_hash.d_1,
                                                            arion_hash.constants_g[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                            arion_hash.constants_h[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                            arion_hash.constants_aff[arion_hash.rounds,:],
                                                            arion_hash.matrix) - next_state_1
        polys_2 = round_function_ArionHash_polynomial_model(current_state_2,
                                                            tmp_2,
                                                            arion_hash.branches,
                                                            arion_hash.d_1,
                                                            arion_hash.constants_g[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                            arion_hash.constants_h[(arion_hash.branches - 1) * (arion_hash.rounds - 1) + 1:(arion_hash.branches - 1) * arion_hash.rounds,:],
                                                            arion_hash.constants_aff[arion_hash.rounds,:],
                                                            arion_hash.matrix) - next_state_2
        push!(polynomials, polys_1[1, 1] - polys_2[1, 1])
        polynomials = append_polynomial_matrix_to_vector(polynomials, polys_1[2:arion_hash.branches, :])
        polynomials = append_polynomial_matrix_to_vector(polynomials, polys_2[2:arion_hash.branches, :])
        if naive_model
            push!(polynomials, tmp_1^arion_hash.e_2 - variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index])
            push!(polynomials, tmp_2^arion_hash.e_2 - variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index])
        else
            push!(polynomials, tmp_1 - variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index]^arion_hash.d_2)
            push!(polynomials, tmp_2 - variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index]^arion_hash.d_2)
        end
    else
        next_state_1 = zero_matrix(P, arion_hash.branches, 1)
        next_state_2 = zero_matrix(P, arion_hash.branches, 1)
        next_state_1[1, 1] = zero(arion_hash.field)
        next_state_2[1, 1] = zero(arion_hash.field)
        for i in 2:arion_hash.branches
            next_state_1[i, 1] = variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
            next_state_2[i, 1] = variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index + i - 1]
        end
        polys_1 = round_function_ArionHash_polynomial_model(current_state_1,
                                                            tmp_1,
                                                            arion_hash.branches,
                                                            arion_hash.d_1,
                                                            arion_hash.constants_g,
                                                            arion_hash.constants_h,
                                                            arion_hash.constants_aff[1,:],
                                                            arion_hash.matrix) - next_state_1
        polys_2 = round_function_ArionHash_polynomial_model(current_state_2,
                                                            tmp_2,
                                                            arion_hash.branches,
                                                            arion_hash.d_1,
                                                            arion_hash.constants_g,
                                                            arion_hash.constants_h,
                                                            arion_hash.constants_aff[1,:],
                                                            arion_hash.matrix) - next_state_2
        push!(polynomials, polys_1[1, 1] - polys_2[1, 1])
        polynomials = append_polynomial_matrix_to_vector(polynomials, polys_1[2:arion_hash.branches, :])
        polynomials = append_polynomial_matrix_to_vector(polynomials, polys_2[2:arion_hash.branches, :])
        if naive_model
            push!(polynomials, tmp_1^arion_hash.e_2 - variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index])
            push!(polynomials, tmp_2^arion_hash.e_2 - variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index])
        else
            push!(polynomials, tmp_1 - variables_1[(arion_hash.rounds - 1) * arion_hash.branches + index]^arion_hash.d_2)
            push!(polynomials, tmp_2 - variables_2[(arion_hash.rounds - 1) * arion_hash.branches + index]^arion_hash.d_2)
        end
    end

    if field_equations
        fes = generate_field_equations(variables)
        polynomials = append_polynomial_matrix_to_vector(polynomials, matrix(fes))
    end
    return polynomials
end
