# Arion SageMath Reference Implementation
Contains a reference implementation of Arion in [SageMath](https://www.sagemath.org/) and the density experiment for the Arion permutation.

### Installation
**Requirements**: All scripts have been developed and tested with SageMath 9.3.

### Examples
**Usage of Arion Reference Implementation**
````
sage: load("Arion.sage")
sage: field = field = GF(1009)
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
````
sage: load("ArionHash.sage")
sage: field = field = GF(1009)
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