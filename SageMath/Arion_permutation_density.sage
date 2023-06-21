"""
    Density Experiment for the Arion permutation.
"""

class ArionDensity:
    
    def __init__(self,
                 field=GF(11),
                 branches=3,
                 rounds=1,
                 d_1=None,
                 d_2=None,
                 constants_g=None,
                 constants_h=None,
                 constants_aff=None,
                 matrix_type=1):
    
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
            self.d_2 = self.d_1
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
        
        if matrix_type == 1:
            self.matrix = matrix.circulant(vector(self.field, range(1, self.branches + 1)))
        elif matrix_type == 2:
            self.matrix = matrix(self.field, self.branches * [list(range(1, self.branches + 1))])
            for i in range(1, self.branches):
                self.matrix[i, i - 1] += 1
                self.matrix[i, i] -= 1
        else:
            raise Exception("Matrix type not implemented.")
        
        print("Arion parameters")
        print("Prime field: " + str(self.field.characteristic()))
        print("Branches: " + str(self.branches))
        print("Rounds: " + str(self.rounds))
        print("Exponent d_1: " + str(self.d_1))
        print("Exponent d_2: " + str(self.d_2))
        print("Matrix:\n" + str(self.matrix))
        print("Constants for the g_i's: " + str(self.constants_g))
        print("Constants for the h_i's: " + str(self.constants_h))
        print("Affine constants: " + str(self.constants_aff))
    
    def reduce_polynomial(self, polynomial, reduction_polynomials):
        for poly in reduction_polynomials:
            polynomial = polynomial % poly
        return polynomial
    
    def reduce_polynomials(self, polynomials, reduction_polynomials):
        for i in range(0, len(polynomials)):
            polynomials[i] = self.reduce_polynomial(polynomials[i], reduction_polynomials)
        return polynomials
    
    def exponentiation_with_reduction(self, polynomial, exponent, reduction_polynomials):
        out = polynomial
        for i in range(1, exponent):
            out *= polynomial
            out = self.reduce_polynomial(out, reduction_polynomials)
        return out
    
    def evaluate_g_and_h_with_reduction(self, x_in, constants_g, constants_h, reduction_polynomials):
        # Constant term
        out_g = constants_g[1]
        out_h = self.field(0)
        # Linear term
        out_g += constants_g[0] * x_in
        out_h += constants_h * x_in
        # Quadratic term
        x_temp = self.exponentiation_with_reduction(x_in, 2, reduction_polynomials)
        out_g += x_temp
        out_h += x_temp
        return out_g, out_h
    
    def GTDS_with_reduction(self, v_in, constants_g, constants_h, reduction_polynomials):
        v_out = copy(v_in)
        v_out[self.branches - 1] = self.exponentiation_with_reduction(v_out[self.branches - 1], self.e_2, reduction_polynomials)
        sigma = v_in[self.branches - 1] + v_out[self.branches - 1]
        for i in range(self.branches - 2, -1, -1):
            v_out[i] = self.exponentiation_with_reduction(v_out[i], self.d_1, reduction_polynomials)
            g_i, h_i = self.evaluate_g_and_h_with_reduction(sigma, constants_g[i], constants_h[i], reduction_polynomials)
            v_out[i] *= g_i
            v_out[i] += h_i
            v_out[i] = self.reduce_polynomial(v_out[i], reduction_polynomials)
            sigma += v_in[i] + v_out[i]
        return v_out
    
    def affine_layer(self, v_in, round_constants):
        return self.matrix * v_in + round_constants
    
    def round_function(self, v_in, constants_g, constants_h, constants_aff, reduction_polynomials):
        v_out = self.GTDS_with_reduction(v_in, constants_g, constants_h, reduction_polynomials)
        v_out = self.affine_layer(v_out, constants_aff)
        return v_out
    
    def permutation(self, v_in, reduction_polynomials):
        v_out = self.matrix * vector(v_in)
        variables = v_out[0].parent().gens()
        for r in range(0, self.rounds):
            print("Round: " + str(r + 1))
            v_out = self.reduce_polynomials(self.round_function(v_out,
                                                                self.constants_g[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                                self.constants_h[r * (self.branches - 1):(r + 1) * (self.branches - 1)],
                                                                vector(self.field, self.constants_aff[r]),
                                                                reduction_polynomials),
                                            reduction_polynomials)
            density = [len(poly.coefficients()) for poly in v_out]
            print("Density of polynomials: " + str(density))
            degrees = [poly.degree() for poly in v_out]
            print("Degree of polynomials: " + str(degrees))
            degrees = [[poly.degree(var) for var in variables] for poly in v_out]
            print("Univariate degrees of polynomials: " + str(degrees))
        return v_out
    
    def density_of_Arion_permutation(self):
        variables = ["x_" + str(i) for i in range(1, self.branches + 1)]
        P = PolynomialRing(self.field, variables, order="degrevlex")
        variables = [P(var) for var in variables]
        field_equations = [var**self.field.characteristic() - var for var in variables]
        print("Computing density of Arion Permutation.")
        print("Maximal possible density: " + str(self.field.characteristic()**self.branches))
        print("Maximal degree: " + str((self.field.characteristic() - 1) * self.branches))
        print("Maximal univariate degree: " + str(self.field.characteristic() - 1))

        polynomials = self.permutation(variables, field_equations)
        return polynomials
