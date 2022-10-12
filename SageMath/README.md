# Arion SageMath Reference Implementation
Contains a reference implementation of Arion in [SageMath](https://www.sagemath.org/) and the density experiment for the Arion permutation.

## Installation
**Requirements**: All scripts have been developed and tested with SageMath 9.3.

## Examples
**Usage of Arion Reference Implementation**
````Python
sage: load("Arion.sage")
sage: field = GF(1009)
sage: branches = 3
sage: rounds = 6
sage: constants_g = [
 [181, 858],
 [127, 665],
 [508, 9],
 [595, 364],
 [350, 572],
 [47, 581],
 [278, 139],
 [630, 832],
 [513, 780],
 [490, 472],
 [918, 465],
 [466, 299]
 ]
sage: constants_h = [17, 805, 73, 795, 3, 666, 514, 42, 314, 900, 466, 719]
sage: constants_aff = [
 [84, 37, 761],
 [530, 421, 892],
 [80, 519, 275],
 [473, 229, 381],
 [21, 239, 781],
 [171, 208, 597]
 ]
sage: arion = Arion(field=field,
                    branches=branches,
                    rounds=rounds,
                    constants_g=constants_g,
                    constants_h=constants_h,
                    constants_aff=constants_aff)
Arion parameters
Prime field: 1009
Branches: 3
Rounds: 6
Exponent d_1: 5
Exponent d_2: 257
Constants for the g_i's: [[181, 858], [127, 665], [508, 9], [595, 364], [350, 572], [47, 581], [278, 139], [630, 832], [513, 780], [490, 472], [918, 465], [466, 299]]
Constants for the h_i's: [17, 805, 73, 795, 3, 666, 514, 42, 314, 900, 466, 719]
Affine constants: [[84, 37, 761], [530, 421, 892], [80, 519, 275], [473, 229, 381], [21, 239, 781], [171, 208, 597]]

sage: plain = branches * [1]
sage: plain
[1, 1, 1]

sage: key = (rounds + 1) * [branches * [2]]
sage: key
[[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]]

sage: cipher = arion.encrypt(plain, key)
sage: cipher
[962, 708, 674]

sage: arion.decrypt(cipher, key) == plain
True
````
**Usage of ArionHash Reference Implementation**
````Python
sage: load("ArionHash.sage")
sage: field = GF(1009)
sage: branches = 3
sage: rounds = 6
sage: capacity = 1
sage: constants_g = [
 [181, 858],
 [127, 665],
 [508, 9],
 [595, 364],
 [350, 572],
 [47, 581],
 [278, 139],
 [630, 832],
 [513, 780],
 [490, 472],
 [918, 465],
 [466, 299]
 ]
sage: constants_h = [17, 805, 73, 795, 3, 666, 514, 42, 314, 900, 466, 719]
sage: constants_aff = [
 [84, 37, 761],
 [530, 421, 892],
 [80, 519, 275],
 [473, 229, 381],
 [21, 239, 781],
 [171, 208, 597]
 ]
sage: arion_hash = ArionHash(field=field,
                       branches=branches,
                       rounds=rounds,
                       capacity=capacity,
                       constants_g=constants_g,
                       constants_h=constants_h,
                       constants_aff=constants_aff)
ArionHash parameters
Prime field: 1009
Branches: 3
Rounds: 6
Sponge capacity: 1
Exponent d_1: 5
Exponent d_2: 257
Constants for the g_i's: [[181, 858], [127, 665], [508, 9], [595, 364], [350, 572], [47, 581], [278, 139], [630, 832], [513, 780], [490, 472], [918, 465], [466, 299]]
Constants for the h_i's: [17, 805, 73, 795, 3, 666, 514, 42, 314, 900, 466, 719]
Affine constants: [[84, 37, 761], [530, 421, 892], [80, 519, 275], [473, 229, 381], [21, 239, 781], [171, 208, 597]]
Initial value: [0]

sage: plain = 3 * [1]
sage: hash_val = arion_hash.hash(plain)
sage: hash_val
748
````
**Usage of Arion density experiment**
````Python
sage: load("Arion_density_experiment.sage")
sage: field = GF(19)
sage: branches = 3
sage: rounds = 3
sage: constants_g = [[5, 3], [14, 9], [6, 10], [18, 16], [1, 6], [14, 15]]
sage: constants_h = [13, 6, 4, 14, 10, 15]
sage: constants_aff = [[11, 2, 1], [1, 0, 7], [9, 15, 6]]
sage: arion_permutation = ArionDensity(field=field,
                                       branches=branches,
                                       rounds=rounds,
                                       constants_g=constants_g,
                                       constants_h=constants_h,
                                       constants_aff=constants_aff)
Arion parameters
Prime field: 19
Branches: 3
Rounds: 3
Exponent d_1: 5
Exponent d_2: 5
Constants for the g_i's: [[5, 3], [14, 9], [6, 10], [18, 16], [1, 6], [14, 15]]
Constants for the h_i's: [13, 6, 4, 14, 10, 15]
Affine constants: [[11, 2, 1], [1, 0, 7], [9, 15, 6]]

sage: polys = arion_permutation.density_of_Arion_permutation()
Computing density of Arion Permutation.
Maximal possible density: 6859
Round: 1
Density of polynomials: [4353, 4362, 4359]
Round: 2
Density of polynomials: [6514, 6508, 6529]
Round: 3
Density of polynomials: [6512, 6465, 6518]
````
**Usage of Arion polynomial model**
````Python
sage: load("Arion.sage")
sage: load("Arion_polynomial_model.sage")
sage: field = GF(11)
sage: branches = 3
sage: rounds = 2
sage: d_1 = 3
sage: d_2 = 7
sage: constants_g = [[3, 10], [6, 10], [6, 7], [5, 7]]
sage: constants_h = [0, 5, 5, 0]
sage: constants_aff = [[3, 9, 8], [9, 2, 7]]
sage: arion = Arion(field=field,
                    branches=branches,
                    rounds=rounds,
                    d_1=d_1,
                    d_2=d_2,
                    constants_g=constants_g,
                    constants_h=constants_h,
                    constants_aff=constants_aff)
sage: plain = branches * [1]
sage: key = (rounds + 1) * [branches * [2]]
sage: cipher = arion.encrypt(plain, key)
sage: cipher
[7, 3, 6]
sage: polys = generate_Arion_polynomials(field=field,
                                         branches=branches,
                                         rounds=rounds,
                                         d_1=d_1,
                                         d_2=d_2,
                                         constants_g=constants_g,
                                         constants_h=constants_h,
                                         constants_aff=constants_aff,
                                         plain=plain,
                                         cipher=cipher,
                                         field_equations=True)
Arion parameters
Prime field: 11
Branches: 3
Rounds: 2
Exponent d_1: 3
Exponent d_2: 7
Constants for the g_i's: [[3, 10], [6, 10], [6, 7], [5, 7]]
Constants for the h_i's: [0, 5, 5, 0]
Affine constants: [[3, 9, 8], [9, 2, 7]]
Plain text: [1, 1, 1]
Cipher text: [7, 3, 6]
Term order: degrevlex

sage: polys[2 * branches + 1] # d_2 = 7
-z_2^7 + x_3__1

sage: gb = ideal(polys).groebner_basis(algorithm="singular:slimgb")
sage: gb
[y_3^2 - 5*y_3 - 5, z_1 - 2, x_1__1 + 5*y_3 - 5, x_2__1 - 3*y_3 - 4, x_3__1 - 5*y_3 - 5, z_2 - 5*y_3 + 1, y_1 - 4*y_3 - 5, y_2 + 3*y_3 + 3]

sage: ideal(gb).variety()
[{y_3: 3, y_2: 10, y_1: 6, z_2: 3, x_3__1: 9, x_2__1: 2, x_1__1: 1, z_1: 2}, 
 {y_3: 2, y_2: 2, y_1: 2, z_2: 9, x_3__1: 4, x_2__1: 10, x_1__1: 6, z_1: 2}]

sage: polys_naive = generate_Arion_polynomials(field=field,
                                               branches=branches,
                                               rounds=rounds,
                                               d_1=d_1,
                                               d_2=d_2,
                                               constants_g=constants_g,
                                               constants_h=constants_h,
                                               constants_aff=constants_aff,
                                               plain=plain,
                                               cipher=cipher,
                                               field_equations=True,
                                               naive_model=True)
sage: polys_naive[2 * branches + 1] # d_2 = 7, p = 11 => e_2 = 3 since 3 * 7 = 21 = 1 mod p - 1
x_3__1^3 - z_2

sage: gb_naive = ideal(polys_naive).groebner_basis(algorithm="singular:slimgb")
sage: gb_naive
[y_3^2 - 5*y_3 - 5, z_1 - 2, x_1__1 + 5*y_3 - 5, x_2__1 - 3*y_3 - 4, x_3__1 - 5*y_3 - 5, z_2 - 5*y_3 + 1, y_1 - 4*y_3 - 5, y_2 + 3*y_3 + 3]
sage: ideal(gb) == ideal(gb_naive)
True
````
