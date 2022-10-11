"""
    Polynomial model of Arion.
"""
load("..\\GitHub\\sca-research\\Arion\\SageMath\\Arion.sage")
#load("Arion.sage")

def generate_variables(branches, rounds):
    variables = []
    variables.append("z_" + str(1))
    for i in range(0, rounds - 1):
        for j in range(0, branches):
            variables.append("x_" + str(j + 1) + "__" + str(i + 1))
        variables.append("z_" + str(i + 2))
    for i in range(0, branches):
        variables.append("y_" + str(i + 1))
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

def key_addition(v_in, key):
    return v_in + key

def round_function(v_in, v_n_std, key, d_1, matrix, constants_g, constants_h, constants_aff):
    v_out = GTDS(v_in, v_n_std, d_1, constants_g, constants_h)
    v_out = affine_layer(v_out, matrix, constants_aff)
    v_out = key_addition(v_out, key)
    return v_out

def generate_Arion_polynomials(field=GF(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), # BLS12
                               branches=3,
                               rounds=6,
                               d_1=None,
                               d_2=None,
                               constants_g=None,
                               constants_h=None,
                               constants_aff=None,
                               plain=None,
                               key=None,
                               cipher=None,
                               termorder="degrevlex",
                               naive_model=False,
                               field_equations=False):
    
    arion = Arion(field=field,
                  branches=branches,
                  rounds=rounds,
                  d_1=d_1,
                  d_2=d_2,
                  constants_g=constants_g,
                  constants_h=constants_h,
                  constants_aff=constants_aff)
    
    if plain is None:
        plain = [field.random_element() for i in range(0, branches)]
    else:
        if len(plain) != branches:
            raise Exception("Plain text length does not match branch number.")
    
    if cipher is None:
        if key is None:
            key = (rounds + 1) * [[field.random_element() for i in range(0, branches)]]
        else:
            if len(key) != branches:
                raise Exception("Key text length does not match branch number.")
            key = (rounds + 1) * [key]
        cipher = arion.encrypt(plain, key)
    else:
        if len(cipher) != branches:
            raise Exception("Cipher text length does not match branch number.")
    
    print("Plain text: " + str(plain))
    if not key is None:
        print("Key: " + str(key[0]))
    print("Cipher text: " + str(plain))
    print("Term order: " + termorder)
    
    variables = generate_variables(branches, rounds)
    P = PolynomialRing(field, variables, order=termorder)
    variables = [P(var) for var in variables]
    polynomials = []
    
    counter = 0
    key = vector(variables[-branches:])
    current_state = arion.matrix * (vector(plain) + key)
    tmp = current_state[branches - 1]
    current_state[branches - 1] = variables[counter]
    counter += 1
    if rounds > 1:
        next_state = vector(variables[counter:counter + branches])
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               key,
                                               arion.d_1,
                                               arion.matrix,
                                               arion.constants_g[0:branches - 1],
                                               arion.constants_h[0:branches - 1],
                                               vector(field, arion.constants_aff[0])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion.d_2)
        for r in range(1, rounds - 1):
            current_state = vector(variables[counter:counter + branches])
            tmp = current_state[branches - 1]
            counter += branches
            current_state[branches - 1] = variables[counter]
            counter += 1
            next_state = vector(variables[counter:counter + branches])
            polynomials.append(list(round_function(current_state,
                                                   tmp,
                                                   key,
                                                   arion.d_1,
                                                   arion.matrix,
                                                   arion.constants_g[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   arion.constants_h[r * (branches - 1):(r + 1) * (branches - 1)],
                                                   vector(field, arion.constants_aff[r])) - next_state))
            if naive_model:
                polynomials.append(tmp**arion.e_2 - current_state[branches - 1])
            else:
                polynomials.append(tmp - current_state[branches - 1]**arion.d_2)
        next_state = vector(cipher)
        current_state = vector(variables[counter:counter + branches])
        tmp = current_state[branches - 1]
        counter += branches
        current_state[branches - 1] = variables[counter]
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               key,
                                               arion.d_1,
                                               arion.matrix,
                                               arion.constants_g[(rounds - 1) * (branches - 1):rounds * (branches - 1)],
                                               arion.constants_h[(rounds - 1) * (branches - 1):rounds * (branches - 1)],
                                               vector(field, arion.constants_aff[rounds - 1])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion.d_2)
    else:
        next_state = vector(cipher)
        polynomials.append(list(round_function(current_state,
                                               tmp,
                                               key,
                                               arion.d_1,
                                               arion.matrix,
                                               arion.constants_g[0:branches - 1],
                                               arion.constants_h[0:branches - 1],
                                               vector(field, arion.constants_aff[0])) - next_state))
        if naive_model:
            polynomials.append(tmp**arion.e_2 - current_state[branches - 1])
        else:
            polynomials.append(tmp - current_state[branches - 1]**arion.d_2)
    
    polynomials = flatten(polynomials)
    
    if field_equations:
        for var in variables:
            polynomials.append(var**field.characteristic() - var)
    
    return polynomials