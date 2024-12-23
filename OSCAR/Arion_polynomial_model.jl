using Oscar
include("Arion.jl")
include("utilities.jl")


function generate_Arion_variables(branches, rounds)
    variables = String[]
    for i in 1:(rounds - 1)
        for j in 1:branches
            push!(variables, "x_" * string(j) * "__" * string(i))
        end
    end
    for i in 1:branches
        push!(variables, "y_" * string(i))
    end
    for i in 1:rounds
        push!(variables, "z_" * string(i))
    end
    return variables
end

function GTDS_Arion_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
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

function round_function_Arion_polynomial_model(v_in, v_n_std, key, branches, d, constants_g, constants_h, constants_aff, mat)
    v_out = GTDS_Arion_polynomial_model(v_in, v_n_std, branches, d, constants_g, constants_h)
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
    if isnothing(plain)
        plain = zero_matrix(arion.field, arion.branches, 1)
        for i in 1:arion.branches
            plain[i, 1] = rand(arion.field)
        end
        if isnothing(cipher)
            print_key = true
            key = zero_matrix(arion.field, arion.branches, 1)
            for i in 1:arion.branches
                key[i, 1] = rand(arion.field)
            end
            cipher = encrypt(plain, key, arion)
        end
    else
        plain = int_vector_to_field_matrix(plain, arion.field)
        if isnothing(cipher)
            print_key = true
            key = zero_matrix(arion.field, arion.branches, 1)
            for i in 1:arion.branches
                key[i, 1] = rand(arion.field)
            end
            cipher = encrypt(plain, key, arion)
        else
            cipher = int_vector_to_field_matrix(cipher, arion.field)
        end
    end

    println("Plain text: ", plain)
    if print_key
        println("Key: ", key)
    end
    println("Cipher text: ", cipher)
    println("Term order: ", termorder)

    if termorder == "degrevlex"
        P, variables = polynomial_ring(arion.field, generate_Arion_variables(arion.branches, arion.rounds), internal_ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = polynomial_ring(arion.field, generate_Arion_variables(arion.branches, arion.rounds), internal_ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end
    
    polynomials = Vector{typeof(variables[1])}()
    key_variables = matrix(variables[arion.branches * (arion.rounds - 1) + 1:arion.branches * arion.rounds])
    index = arion.branches * arion.rounds + 1
    if arion.rounds > 1
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[index]
        next_state = matrix(variables[1:arion.branches])
        polynomials = round_function_Arion_polynomial_model(current_state,
                                                            tmp,
                                                            key_variables,
                                                            arion.branches,
                                                            arion.d_1,
                                                            arion.constants_g[1:arion.branches,:],
                                                            arion.constants_h[1:arion.branches,:],
                                                            arion.constants_aff[1,:],
                                                            arion.matrix) - next_state
        polynomials = [polynomials; vec(polys[:, 1])]
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion.d_2)
        end
        index += 1
        for r in 1:arion.rounds - 2
            current_state = matrix(variables[(r - 1) * arion.branches + 1:r * arion.branches])
            tmp = current_state[arion.branches, 1]
            current_state[arion.branches, 1] = variables[index]
            next_state = matrix(variables[r * arion.branches + 1:(r + 1) * arion.branches])
            polynomials = round_function_Arion_polynomial_model(current_state,
                                                                tmp,
                                                                key_variables,
                                                                arion.branches,
                                                                arion.d_1,
                                                                arion.constants_g[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                arion.constants_h[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                arion.constants_aff[r + 1,:],
                                                                arion.matrix) - next_state
            polynomials = [polynomials; vec(polys[:, 1])]
            if naive_model
                push!(polynomials, tmp^arion.e_2 - variables[index])
            else
                push!(polynomials, tmp - variables[index]^arion.d_2)
            end
            index += 1
        end
        current_state = matrix(variables[(arion.rounds - 2) * arion.branches + 1:(arion.rounds - 1) * arion.branches])
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[index]
        polynomials = round_function_Arion_polynomial_model(current_state,
                                                            tmp,
                                                            key_variables,
                                                            arion.branches,
                                                            arion.d_1,
                                                            arion.constants_g[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                            arion.constants_h[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                            arion.constants_aff[arion.rounds,:],
                                                            arion.matrix) - cipher
        polynomials = [polynomials; vec(polys[:, 1])]
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion.d_2)
        end
    else
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[1]
        polynomials = round_function_Arion_polynomial_model(current_state,
                                                            tmp,
                                                            key_variables,
                                                            arion.branches,
                                                            arion.d_1,
                                                            arion.constants_g,
                                                            arion.constants_h,
                                                            arion.constants_aff[1,:],
                                                            arion.matrix) - cipher
        polynomials = [polynomials; vec(polys[:, 1])]
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion.d_2)
        end
    end
    
    if field_equations
        fes = generate_field_equations(variables)
        polynomials = [polynomials; fes]
    end
    return polynomials
end
