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

function generate_input_output_variables(branches)
    input_variables = String[]
    output_variables = String[]
    for i in 1:branches
        push!(input_variables, "x_in__" * string(i))
        push!(output_variables, "x_out__" * string(i))
    end
    return input_variables, output_variables
end

function GTDS_Arion_polynomial_model(v_in, 
                                     v_n_std, 
                                     branches, 
                                     d, 
                                     constants_g, 
                                     constants_h;
                                     sigma_type="I")
    v_out = zero_matrix(base_ring(v_in), branches, 1)
    v_out[branches, 1] = v_in[branches, 1]
        sigma = v_n_std + v_out[branches, 1]
    for i in branches - 1:-1:1
        v_out[i, 1] = v_in[i, 1]^d
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i,:], constants_h[i,:])
        v_out[i, 1] *= g_i
        v_out[i, 1] += h_i
        if sigma_type == "I"
            sigma += v_in[i, 1] + v_out[i, 1]
        elseif sigma_type == "II"
            sigma += sigma^2 + v_in[i, 1]
        elseif sigma_type == "III"
            sigma += v_in[i, 1]
        else
            println("Given sigma Type not implemented!")
            return
        end
    end
    return v_out
end

function round_function_Arion_polynomial_model(v_in, 
                                               v_n_std, 
                                               key, 
                                               branches, 
                                               d, 
                                               constants_g, 
                                               constants_h, 
                                               constants_aff, 
                                               mat;
                                               sigma_type="I")
    v_out = GTDS_Arion_polynomial_model(v_in, 
                                        v_n_std, 
                                        branches, 
                                        d, 
                                        constants_g, 
                                        constants_h, 
                                        sigma_type=sigma_type)
    v_out = affine_layer(v_out, constants_aff, mat)
    v_out = key_addition(v_out, key)
    return v_out
end

function generate_generic_Arion_polynomials(;arion=Arion_constructor(),
                                             sigma_type="I",
                                             termorder="degrevlex",
                                             naive_model=false,
                                             info_level=0)
    if info_level > 0
        println("Sigma type in Arion GTDS: ", sigma_type)
        println("Term order: ", termorder)
    end

    variables = generate_Arion_variables(arion.branches, arion.rounds)
    input_variables, output_variables = generate_input_output_variables(arion.branches)
    variables = [variables; input_variables; output_variables]
    if termorder == "degrevlex"
        P, variables = PolynomialRing(arion.field, variables, ordering=:degrevlex)
    elseif termorder == "lex"
        P, variables = PolynomialRing(arion.field, variables, ordering=:lex)
    else
        println("Term order ", termorder, " is not implemented.")
        return
    end

    l = length(variables)
    plain = matrix(variables[l - 2 * arion.branches + 1:l - arion.branches])
    cipher = matrix(variables[l - arion.branches + 1:l])
    variables = variables[1:l - 2 * arion.branches]

    polynomials = Vector{typeof(variables[1])}()
    key_variables = matrix(variables[arion.branches * (arion.rounds - 1) + 1:arion.branches * arion.rounds])
    if arion.rounds > 1
        index = arion.branches * arion.rounds + 1
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[index]
        next_state = matrix(variables[1:arion.branches])
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_Arion_polynomial_model(current_state,
                                                                                               tmp,
                                                                                               key_variables,
                                                                                               arion.branches,
                                                                                               arion.d_1,
                                                                                               arion.constants_g[1:arion.branches,:],
                                                                                               arion.constants_h[1:arion.branches,:],
                                                                                               arion.constants_aff[1,:],
                                                                                               arion.matrix,
                                                                                               sigma_type=sigma_type) - next_state)
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
            polynomials = append_polynomial_matrix_to_vector(polynomials,
                                                             round_function_Arion_polynomial_model(current_state,
                                                                                                   tmp,
                                                                                                   key_variables,
                                                                                                   arion.branches,
                                                                                                   arion.d_1,
                                                                                                   arion.constants_g[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                                                   arion.constants_h[(arion.branches - 1) * r + 1:(arion.branches - 1) * (r + 1),:],
                                                                                                   arion.constants_aff[r + 1,:],
                                                                                                   arion.matrix,
                                                                                                   sigma_type=sigma_type) - next_state)
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
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_Arion_polynomial_model(current_state,
                                                                                               tmp,
                                                                                               key_variables,
                                                                                               arion.branches,
                                                                                               arion.d_1,
                                                                                               arion.constants_g[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                                                               arion.constants_h[(arion.branches - 1) * (arion.rounds - 1) + 1:(arion.branches - 1) * arion.rounds,:],
                                                                                               arion.constants_aff[arion.rounds,:],
                                                                                               arion.matrix,
                                                                                               sigma_type=sigma_type) - cipher)
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[index])
        else
            push!(polynomials, tmp - variables[index]^arion.d_2)
        end
    else
        current_state = arion.matrix * (plain + key_variables)
        tmp = current_state[arion.branches, 1]
        current_state[arion.branches, 1] = variables[arion.branches + 1]
        polynomials = append_polynomial_matrix_to_vector(polynomials, 
                                                         round_function_Arion_polynomial_model(current_state,
                                                                                               tmp,
                                                                                               key_variables,
                                                                                               arion.branches,
                                                                                               arion.d_1,
                                                                                               arion.constants_g,
                                                                                               arion.constants_h,
                                                                                               arion.constants_aff[1,:],
                                                                                               arion.matrix,
                                                                                               sigma_type=sigma_type) - cipher)
        if naive_model
            push!(polynomials, tmp^arion.e_2 - variables[arion.branches + 1])
        else
            push!(polynomials, tmp - variables[arion.branches + 1]^arion.d_2)
        end
    end
    
    return polynomials
end
