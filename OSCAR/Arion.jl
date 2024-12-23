using Oscar
include("utilities.jl")

struct Arion
    field::FqField
    branches::Int64
    rounds::Int64
    d_1::Int64
    e_1::Int64
    d_2::Int64
    e_2::Int64
    matrix::FqMatrix
    constants_g::FqMatrix
    constants_h::FqMatrix
    constants_aff::FqMatrix
end

function Arion_constructor(;field=GF(1009),
                            branches=3,
                            rounds=3,
                            d_1=-1,
                            d_2=-1,
                            constants_g=nothing,
                            constants_h=nothing,
                            constants_aff=nothing,
                            matrix_type=1)
            
    q = order(field)
    
    if d_1 == -1
        counter = 2
        while gcd(counter, q - 1) != 1
            counter += 1
        end
        d_1 = counter
    else
        if gcd(d_1, q - 1) != 1
            print(d_1, " does not induce permutation over ", field)
            return
        end
    end
    e_1 = gcdx(d_1, q - 1)[2]
    while e_1 < 0
        e_1 += q - 1
    end
    while q - 1 < e_1 
        e_1 -= q - 1
    end
    if d_2 == -1
        found = false
        for d in [257, 161, 129, 125, 123, 121]
            if gcd(d, q - 1) == 1
                d_2 = d
                found = true
                break
            end
        end
        if !found
            println("Could not find power permutation d_2 with exponent in ", [257, 161, 129, 125, 123, 121],  + " over " + field)
        end
    else
        if gcd(d_2, q - 1) != 1
            print(d_2, " does not induce permutation over ", field)
            return
        end
    end
    e_2 = gcdx(d_2, q - 1)[2]
    while e_2 < 0
        e_2 += q - 1
    end
    while q - 1 < e_2 
        e_2 -= q - 1
    end

    mat = zero_matrix(field, branches, branches)
    if matrix_type == 1
        row = circshift(Vector(1:branches), -1)
        for i in 1:branches
            row = circshift(row, 1)
            for j in 1:branches
                mat[i, j] = field(row[j])
            end
    end
    elseif matrix_type == 2
        current_row = Vector(1:branches)
        for i in 1:branches
            for j in 1:branches
                mat[i, j] = field(current_row[j])
            end
            if i > 1
                mat[i, i - 1] += field(1)
                mat[i, i] -= field(1)
            end
        end
    else
        println("Matrix type not implemented.")
        return
    end

    if isnothing(constants_g)
        constants_g = zero_matrix(field, (branches - 1) * rounds, 2)
        counter = 1
        for i in 1:rounds
            for j in 1:(branches - 1)
                c_1 = rand(field)
                c_2 = rand(field)
                while legendre_symbol(c_1^2 - 4 * c_2, q) != -1
                    c_1 = rand(field)
                    c_2 = rand(field)
                end
                constants_g[counter, 1] = c_1
                constants_g[counter, 2] = c_2
                counter += 1
            end
        end
    else
        constants_g = matrix(field, reshape(constants_g, (branches - 1) * rounds, 2))
        if ncols(constants_g) * nrows(constants_g) != rounds * (branches - 1) * 2
            println("Number of round constants for the polynomials g_i does not match branches minus one times the number of rounds times two.")
            return
        end
    end

    if isnothing(constants_h)
        constants_h = zero_matrix(field, (branches - 1) * rounds, 1)
        counter = 1
        for i in 1:rounds
            for j in 1:(branches - 1)
                constants_h[counter, 1] = rand(field)
                counter += 1
            end
        end
    else
        constants_h = matrix(field, reshape(constants_h, (branches - 1) * rounds, 1))
        if ncols(constants_h) * nrows(constants_h) != rounds * (branches - 1)
            println("Number of round constants for the polynomials h_i does not match branches minus one times the number of rounds times two.")
            return
        end
    end

    if isnothing(constants_aff)
        constants_aff = zero_matrix(field, rounds, branches)
        for i in 1:rounds
            for j in 1:branches
                constants_aff[i, j] = rand(field)
            end
        end
    else
        constants_aff = matrix(field, reshape(constants_aff, rounds, branches))
        if ncols(constants_aff) * nrows(constants_aff) != rounds * branches
            println("Number of affine round constants does not match branches times the number of rounds.")
            return
        end
    end

    println("Arion parameters")
    println("Prime field: ", q)
    println("Branches: ", branches)
    println("Rounds: ", rounds)
    println("Exponent d_1: ", d_1)
    println("Exponent d_2: ", d_2)
    println("Matrix: ", mat)
    println("Constants for the g_i's: ", constants_g)
    println("Constants for the h_i's: ", constants_h)
    println("Affine constants: ", constants_aff)

    return Arion(field, branches, rounds, d_1, e_1, d_2, e_2, mat, constants_g, constants_h, constants_aff)
end

function evaluate_g_and_h(x_in, constants_g, constants_h)
    # Constant term
    out_g = constants_g[2]
    out_h = 0
    # Linear term
    x_temp = x_in
    out_g += constants_g[1, 1] * x_temp
    out_h += constants_h[1, 1] * x_temp
    # Quadratic term
    x_temp *= x_in
    out_g += x_temp
    out_h += x_temp
    return out_g, out_h
end

function GTDS(v_in, branches, d, e, constants_g, constants_h)
    v_out = zero_matrix(base_ring(v_in), branches, 1)
    v_out[branches, 1] = v_in[branches, 1]^e
    sigma = v_in[branches, 1] + v_out[branches, 1]
    for i in branches - 1:-1:1
        v_out[i, 1] = v_in[i, 1]^d
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i,:], constants_h[i,:])
        v_out[i, 1] *= g_i
        v_out[i, 1] += h_i
        sigma += v_in[i, 1] + v_out[i, 1]
    end
    return v_out
end

function affine_layer(v_in, constants_aff, mat)
    v_out = mat * v_in
    return v_out + transpose(constants_aff)
end

function key_addition(plain, key)
    return plain + key
end

function round_function(v_in, key, branches, d, e, constants_g, constants_h, constants_aff, mat)
    v_out = GTDS(v_in, branches, d, e, constants_g, constants_h)
    v_out = affine_layer(v_out, constants_aff, mat)
    v_out = key_addition(v_out, key)
    return v_out
end

function encrypt(plain, key, arion::Arion)
    cipher = arion.matrix * key_addition(plain, key)
    for r in 1:arion.rounds
        cipher = round_function(cipher,
                                key,
                                arion.branches,
                                arion.d_1,
                                arion.e_2,
                                arion.constants_g[(r - 1) * (arion.branches - 1) + 1:r * (arion.branches - 1), :],
                                arion.constants_h[(r - 1) * (arion.branches - 1) + 1:r * (arion.branches - 1), :],
                                arion.constants_aff[r, :],
                                arion.matrix)
    end
    return cipher
end
