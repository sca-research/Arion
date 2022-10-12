"""
    Polynomial model of arion_hash.
"""

load("ArionHash.sage")

def generate_variables(branches, rate, rounds):
    variables = []
    for i in range(0, rate):
        variables.append("x_in_" + str(i + 1))
    variables.append("z_" + str(1))
    for i in range(0, rounds - 1):
        for j in range(0, branches):
            variables.append("x_" + str(j + 1) + "__" + str(i + 1))
        variables.append("z_" + str(i + 2))
    for i in range(1, branches):
        variables.append("x_out_" + str(i + 1))
    return variables

def evaluate_g_and_h(x_in, constants_g, constant_h):
    # Constant term
    out_g = constants_g[1]
    out_h = field(0)
    # Linear term
    out_g += x_in * constants_g[0]
    out_h += x_in * constant_h
    # Quadratic term
    x_temp = x_in**2
    out_g += x_temp
    out_h += x_temp
    return out_g, out_h

def GTDS(v_in, v_n_std, d_1, constants_g, constants_h):
    branches = len(v_in)
    v_out = copy(v_in)
    sigma = v_n_std + v_out[branches - 1]
    for i in range(branches - 2, -1, -1):
        v_out[i] = v_out[i]**d_1
        g_i, h_i = evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
        v_out[i] *= g_i
        v_out[i] += h_i
        sigma += v_in[i] + v_out[i]
    return v_out

def affine_layer(v_in, matrix, round_constants):
    return matrix * v_in + round_constants

def round_function(v_in, v_n_std, d_1, matrix, constants_g, constants_h, constants_aff):
    v_out = GTDS(v_in, v_n_std, d_1, constants_g, constants_h)
    v_out = affine_layer(v_out, matrix, constants_aff)
    return v_out

def generate_ArionHash_polynomials(field=GF(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), # BLS12
                                   branches=3,
                                   rounds=6,
                                   capacity=2,
                                   d_1=None,
                                   d_2=None,
                                   constants_g=None,
                                   constants_h=None,
                                   constants_aff=None,
                                   hash_val=None,
                                   termorder="degrevlex",
                                   naive_model=False,
                                   field_equations=False):
    
    arion_hash = ArionHash(field=field,
                           branches=branches,
                           rounds=rounds,
                           capacity=capacity,
                           d_1=d_1,
                           d_2=d_2,
                           constants_g=constants_g,
                           constants_h=constants_h,
                           constants_aff=constants_aff)
    print_plain = False
    if hash_val is None:
        plain = field.random_element()
        hash_val = arion_hash.hash(plain)
        print_plain = True
    
    if print_plain:
        print("Plain text: " + str(plain))
    print("Hash value: " + str(hash_val))
    print("Term order: " + termorder)
    
    variables = generate_variables(branches, arion_hash.rate, rounds)
    P = PolynomialRing(field, variables, order=termorder)
    variables = [P(var) for var in variables]
    polynomials = []
    counter = arion_hash.rate
    
    current_state = vector(variables[0:counter] + arion_hash.initial_value)
    current_state = arion_hash.matrix * current_state
    tmp = current_state[branches - 1]
    current_state[branches - 1] = variables[counter]
    counter += 1
    if rounds > 1:
        next_state = vector(variables[counter:counter + branches])
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               arion_hash.d_1,
                                               arion_hash.matrix,
                                               arion_hash.constants_g[0:branches - 1],
                                               arion_hash.constants_h[0:branches - 1],
                                               vector(field, arion_hash.constants_aff[0])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion_hash.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion_hash.d_2)
        for r in range(1, rounds - 1):
            current_state = vector(variables[counter:counter + branches])
            tmp = current_state[branches - 1]
            counter += branches
            current_state[branches - 1] = variables[counter]
            counter += 1
            next_state = vector(variables[counter:counter + branches])
            polynomials.append(list(round_function(current_state,
                                                   tmp,
                                                   arion_hash.d_1,
                                                   arion_hash.matrix,
                                                   arion_hash.constants_g[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   arion_hash.constants_h[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   vector(field, arion_hash.constants_aff[r])) - next_state))
            if naive_model:
                polynomials.append(tmp**arion_hash.e_2 - current_state[branches - 1])
            else:
                polynomials.append(tmp - current_state[branches - 1]**arion_hash.d_2)
        current_state = vector(variables[counter:counter + branches])
        tmp = current_state[branches - 1]
        counter += branches
        current_state[branches - 1] = variables[counter]
        next_state = vector([hash_val] + variables[-branches + 1:])
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               arion_hash.d_1,
                                               arion_hash.matrix,
                                               arion_hash.constants_g[(rounds - 1) * (branches - 1):rounds * (branches - 1)],
                                               arion_hash.constants_h[(rounds - 1) * (branches - 1):rounds * (branches - 1)],
                                               vector(field, arion_hash.constants_aff[rounds - 1])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion_hash.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion_hash.d_2)
    else:
        next_state = vector([hash_val] + variables[-branches + 1:])
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               arion_hash.d_1,
                                               arion_hash.matrix,
                                               arion_hash.constants_g[0:branches - 1],
                                               arion_hash.constants_h[0:branches - 1],
                                               vector(field, arion_hash.constants_aff[0])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion_hash.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion_hash.d_2)
    
    polynomials = flatten(polynomials)
    
    if field_equations:
        for var in variables:
            polynomials.append(var**field.characteristic() - var)
    
    return polynomials
