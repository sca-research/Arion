"""
    Reference Implementatopm of Arion.
"""

class Arion:
    
    def __init__(self,
                 field=GF(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001), # BLS12
                 branches=3,
                 rounds=6,
                 d_1=None,
                 d_2=None,
                 constants_g=None,
                 constants_h=None,
                 constants_aff=None):
        self.field = field
        if not self.field.is_prime_field():
            raise Exception("Only prime fields are implemented. " + str(self.field) + " is not a prime field.")
        self.branches = branches
        self.rounds = rounds
        
        self.d_1 = d_1
        if self.d_1 is None:
            self.d_1 = 2
            while gcd(self.d_1, self.field.characteristic() - 1) != 1:
                self.d_1 += 1
        else:
            if gcd(self.d_1, self.field.characteristic() - 1) != 1:
                raise Exception("d_1 = " + str(self.d_1) + " does not induce a power permutation over " + str(self.field) + ".")
        self.e_1 = xgcd(self.d_1, self.field.characteristic() - 1)[1] % (self.field.characteristic() - 1)
        
        self.d_2 = d_2
        if self.d_2 is None:
            found = False
            for d in [257, 161, 129, 125, 123, 121]:
                if gcd(d, self.field.characteristic() - 1) == 1:
                    self.d_2 = d
                    found = True
                    break
            if not found:
                raise Exception("Could not find power permutation d_2 with exponent in " + str([257, 161, 129, 125, 123, 121])  + " over " + str(self.field) + ".")
        else:
            if gcd(self.d_2, self.field.characteristic() - 1) != 1:
                raise Exception("d_2 = " + str(self.d_2) + " does not induce a power permutation over " + str(self.field) + ".")
        self.e_2 = xgcd(self.d_2, self.field.characteristic() - 1)[1] % (self.field.characteristic() - 1)
        
        self.constants_g = constants_g
        if self.constants_g is None:
            self.constants_g = []
            for r in range(0, self.rounds):
                for k in range(0, self.branches - 1):
                    c_1 = self.field.random_element()
                    c_2 = self.field.random_element()
                    while legendre_symbol(c_1**2 - 4 * c_2, self.field.characteristic()) != -1:
                        c_1 = self.field.random_element()
                        c_2 = self.field.random_element()
                    self.constants_g.append([c_1, c_2])
        else:
            if len(flatten(self.constants_g)) != self.rounds * (self.branches - 1) * 2:
                raise Exception("Number of round constants for the polynomials g_i does not match branches minus one times the number of rounds times two.")
            for pair in self.constants_g:
                if legendre_symbol(pair[0]**2 - 4 * pair[1], self.field.characteristic()) != -1:
                    raise Exception("The constant pair (" + str(pair[0]) + ", " + str(pair[1]) + ") for g_i has zeros over " + str(self.field) + ".")
        
        self.constants_h = constants_h
        if self.constants_h is None:
            self.constants_h = []
            for i in range(0, self.rounds * (self.branches - 1)):
                self.constants_h.append(self.field.random_element())
        else:
            if len(self.constants_h) != self.rounds * (self.branches - 1):
                raise Exception("Number of round constants for the polynomials h_i does not match branches minus one times the number of rounds.")
        
        self.constants_aff = constants_aff
        if self.constants_aff is None:
            self.constants_aff = []
            for r in range(0, self.rounds):
                self.constants_aff.append([self.field.random_element() for i in range(0, self.branches)])
        else:
            if len(flatten(self.constants_aff))!= self.branches * self.rounds:
                raise Exception("Number of affine round constants does not match branches times the number of rounds.")
        
        current_row = list(range(1, self.branches + 1))
        self.matrix = [current_row]
        for i in range(1, self.branches):
            current_row = current_row[-1:] + current_row[:-1]
            self.matrix.append(current_row)
        self.matrix = matrix(self.field, self.matrix)
        try:
            self.matrix_inv = self.matrix.inverse()
        except:
            raise Exception("Circulant matrix circ(1, ..., " + str(self.branches) + ") is not invertible over " + str(self.field) + ".")
        
        print("Arion parameters")
        print("Prime field: " + str(self.field.characteristic()))
        print("Branches: " + str(self.branches))
        print("Rounds: " + str(self.rounds))
        print("Exponent d_1: " + str(self.d_1))
        print("Exponent d_2: " + str(self.d_2))
        print("Constants for the g_i's: " + str(self.constants_g))
        print("Constants for the h_i's: " + str(self.constants_h))
        print("Affine constants: " + str(self.constants_aff))
    
    def evaluate_g_and_h(self, x_in, constants_g, constant_h):
        # Constant term
        out_g = constants_g[1]
        out_h = self.field(0)
        # Linear term
        out_g += x_in * constants_g[0]
        out_h += x_in * constant_h
        # Quadratic term
        x_temp = x_in**2
        out_g += x_temp
        out_h += x_temp
        return out_g, out_h
    
    def GTDS(self, v_in, constants_g, constants_h):
        v_out = copy(v_in)
        v_out[self.branches - 1] = v_out[self.branches - 1]**self.e_2
        sigma = v_in[self.branches - 1] + v_out[self.branches - 1]
        for i in range(self.branches - 2, -1, -1):
            v_out[i] = v_out[i]**self.d_1
            g_i, h_i = self.evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
            v_out[i] *= g_i
            v_out[i] += h_i
            sigma += v_in[i] + v_out[i]
        return v_out
    
    def affine_layer(self, v_in, round_constants):
        return self.matrix * v_in + round_constants
    
    def key_addition(self, v_in, key):
        return v_in + key
    
    def encryption_round_function(self, v_in, key, constants_g, constants_h, constants_aff):
        v_out = self.GTDS(v_in, constants_g, constants_h)
        v_out = self.affine_layer(v_out, constants_aff)
        v_out = self.key_addition(v_out, key)
        return v_out
    
    def encrypt(self, plain, key):
        if len(flatten(key)) != (self.rounds + 1) * self.branches:
            raise Exception("Number of keys does not matach rounds plus one times branches.")
        
        cipher = self.key_addition(vector(self.field, plain), vector(self.field, key[0]))
        cipher = self.matrix * cipher
        for r in range(0, self.rounds):
            cipher = self.encryption_round_function(cipher,
                                                    vector(self.field, key[r + 1]),
                                                    self.constants_g[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                    self.constants_h[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                    vector(self.field, self.constants_aff[r]))
        return list(cipher)
    
    def GTDS_inverse(self, v_in, constants_g, constants_h):
        v_out = copy(v_in)
        v_out[self.branches - 1] = v_out[self.branches - 1]**self.d_2
        sigma = v_in[self.branches - 1] + v_out[self.branches - 1]
        for i in range(self.branches - 2, -1, -1):
            g_i, h_i = self.evaluate_g_and_h(sigma, constants_g[i], constants_h[i])
            v_out[i] -= h_i
            v_out[i] *= g_i**(self.field.characteristic() - 2)
            v_out[i] = v_out[i]**self.e_1
            sigma += v_in[i] + v_out[i]
        return v_out
    
    def inverse_affine_layer(self, v_in, round_constants):
        return self.matrix_inv * (v_in - round_constants)
    
    def key_substraction(self, v_in, key):
        return v_in - key
    
    def decryption_round_function(self, v_in, key, constants_g, constants_h, constants_aff):
        v_out = self.key_substraction(v_in, key)
        v_out = self.inverse_affine_layer(v_out, constants_aff)
        v_out = self.GTDS_inverse(v_out, constants_g, constants_h)
        return v_out
    
    def decrypt(self, plain, key):
        if len(flatten(key)) != (self.rounds + 1) * self.branches:
            raise Exception("Number of keys does not matach rounds plus one times branches.")
        
        v_out = vector(self.field, plain)
        for r in range(self.rounds - 1, -1, -1):
            v_out = self.decryption_round_function(v_out,
                                                   vector(self.field, key[r + 1]),
                                                   self.constants_g[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                   self.constants_h[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                   vector(self.field, self.constants_aff[r]))
        v_out = self.matrix_inv * v_out
        v_out = self.key_substraction(v_out, vector(self.field, key[0]))
        return list(v_out)
