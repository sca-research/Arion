# Arion Gröbner Basis Analysis
Contains the [OSCAR](https://oscar.computeralgebra.de/) Gröbner basis attack implementation of Arion and ArionHash.

## Installation
**Requirements**: All scripts have been developed and tested with [Julia](https://julialang.org/) 1.8.4 and [OSCAR](https://oscar.computeralgebra.de/) 0.12.0.

Currently, [OSCAR](https://oscar.computeralgebra.de/) can only be installed on Linux systems.
Therefore, Windows users first have to activate the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install).

To install [Julia](https://julialang.org/) one can for example build it from source via
```shell
git clone https://github.com/JuliaLang/julia.git
cd julia
git checkout v1.8.4
make
./julia
```
To install [OSCAR](https://oscar.computeralgebra.de/) in [Julia](https://julialang.org/) execute
```julia
julia> using Pkg
julia> Pkg.add(name="Oscar", version="0.12.0")
julia> using Oscar
```
For parallelization the following [Julia](https://julialang.org/) must also be installed
```julia
julia> using Pkg
julia> Pkg.add("Folds")
```


## Examples
**Usage of Arion polynomial model in OSCAR**
```julia
julia> include("Arion_polynomial_model.jl")
julia> field = GF(17)
julia> branches = 3
julia> rounds = 2
julia> d_1 = 3
julia> d_2 = 5
julia> constants_g = [14 12; 16 7; 15 4; 7 2]
julia> constants_h = [3; 9; 4; 3]
julia> constants_aff = [14 12 8; 0 12 13]
julia> arion = Arion_constructor(field=field,
                                 branches=branches,
                                 rounds=rounds,
                                 d_1=d_1,
                                 d_2=d_2,
                                 constants_g=constants_g,
                                 constants_h=constants_h,
                                 constants_aff=constants_aff);
julia> plain = [13; 0; 14]
julia> key = [12; 7; 12]
julia> cipher = encrypt(plain, key, arion)
[1; 5; 3]
julia> polys = generate_Arion_polynomials(arion=arion,
                                          plain=plain,
                                          cipher=cipher,
                                          field_equations=true);
julia> gb = groebner_basis_f4(ideal(polys), nr_thrds=16, info_level=2)
Gröbner basis with elements
1 -> z_2 + 15
2 -> z_1 + 3
3 -> y_3 + 5
4 -> y_2 + 10
5 -> y_1 + 5
6 -> x_3__1 + 2
7 -> x_2__1 + 11
8 -> x_1__1 + 15
with respect to the ordering 
degrevlex([x_1__1, x_2__1, x_3__1, y_1, y_2, y_3, z_1, z_2])

julia> polys_naive = generate_Arion_polynomials(arion=arion,
                                                plain=plain,
                                                cipher=cipher,
                                                field_equations=true,
                                                naive_model=true);
julia> gb_naive = groebner_basis_f4(ideal(polys_naive), nr_thrds=16, info_level=2)
Gröbner basis with elements
1 -> z_2 + 15
2 -> z_1 + 3
3 -> y_3 + 5
4 -> y_2 + 10
5 -> y_1 + 5
6 -> x_3__1 + 2
7 -> x_2__1 + 11
8 -> x_1__1 + 15
with respect to the ordering 
degrevlex([x_1__1, x_2__1, x_3__1, y_1, y_2, y_3, z_1, z_2])
```

**Usage of ArionHash polynomial model in OSCAR**
```julia
julia> include("ArionHash_polynomial_model.jl")
julia> field = GF(17)
julia> branches = 3
julia> rounds = 2
julia> capacity = 2
julia> d_1 = 3
julia> d_2 = 5
julia> constants_g = [14 12; 16 7; 15 4; 7 2]
julia> constants_h = [3; 9; 4; 3]
julia> constants_aff = [14 12 8; 0 12 13]
julia> arion_hash = ArionHash_constructor(field=field,
                                          branches=branches,
                                          rounds=rounds,
                                          capacity=capacity,
                                          d_1=d_1,
                                          d_2=d_2,
                                          constants_g=constants_g,
                                          constants_h=constants_h,
                                          constants_aff=constants_aff);
julia> plain = [5]
julia> hash_val = hash(plain, arion_hash)
7
julia> polys = generate_ArionHash_polynomials(arion_hash=arion_hash,
                                              hash_val=hash_val,
                                              field_equations=true);
julia> gb = groebner_basis_f4(ideal(polys), nr_thrds=16, info_level=2)
Gröbner basis with elements
1 -> z_2 + 16
2 -> z_1 + 6
3 -> x_out__3 + 5
4 -> x_out__2 + 7
5 -> x_3__1 + 16
6 -> x_2__1
7 -> x_1__1 + 12
8 -> x_in__1 + 12
with respect to the ordering
degrevlex([x_in__1, x_1__1, x_2__1, x_3__1, x_out__2, x_out__3, z_1, z_2])

julia> polys_naive = generate_ArionHash_polynomials(arion_hash=arion_hash,
                                                    hash_val=hash_val,
                                                    field_equations=true,
                                                    naive_model=true);
julia> gb_naive = groebner_basis_f4(ideal(polys_naive), nr_thrds=16, info_level=2)
Gröbner basis with elements
1 -> z_2 + 16
2 -> z_1 + 6
3 -> x_out__3 + 5
4 -> x_out__2 + 7
5 -> x_3__1 + 16
6 -> x_2__1
7 -> x_1__1 + 12
8 -> x_in__1 + 12
with respect to the ordering
degrevlex([x_in__1, x_1__1, x_2__1, x_3__1, x_out__2, x_out__3, z_1, z_2])
```

**Usage of ArionHash collision polynomial model in OSCAR**
```julia
julia> include("ArionHash_polynomial_model.jl")
julia> field = GF(11)
julia> branches = 3
julia> rounds = 2
julia> capacity = 2
julia> d_1 = 3
julia> d_2 = 3
julia> constants_g = [3 1; 5 7; 4 4; 7 2]
julia> constants_h = [3; 9; 4; 3]
julia> constants_aff = [3 1 8; 0 1 2]
julia> arion_hash = ArionHash_constructor(field=field,
                                          branches=branches,
                                          rounds=rounds,
                                          capacity=capacity,
                                          d_1=d_1,
                                          d_2=d_2,
                                          constants_g=constants_g,
                                          constants_h=constants_h,
                                          constants_aff=constants_aff);
julia> polys = generate_ArionHash_collision_polynomials(arion_hash=arion_hash,
                                                        field_equations=true);
julia> variables = gens(parent(polys[1]));
julia> I = ideal(polys);
julia> J = ideal([variables[length(variables)] - 4]) # guess of one output state variable
julia> gb = groebner_basis_f4(I + J, nr_thrds=16, info_level=2)
Gröbner basis with elements
1 -> z_2_2 + 7
2 -> z_2_1 + 5
3 -> x_out_2__3 + 3
4 -> x_out_2__2 + 10
5 -> x_2_3__1 + 2
6 -> x_2_2__1 + 8
7 -> x_2_1__1
8 -> x_in_2__1 + 2
9 -> z_1_1 + z_1_2 + 1
10 -> x_out_1__3 + 7*z_1_2 + 8
11 -> x_out_1__2 + 6*z_1_2 + 8
12 -> x_1_3__1 + 6*z_1_2
13 -> x_1_2__1 + 9*z_1_2 + 5
14 -> x_1_1__1 + 4*z_1_2 + 6
15 -> x_in_1__1 + 10*z_1_2 + 6
16 -> z_1_2^2 + 7*z_1_2
with respect to the ordering
degrevlex([x_in_1__1, x_1_1__1, x_1_2__1, x_1_3__1, x_out_1__2, x_out_1__3, z_1_1, z_1_2, x_in_2__1, x_2_1__1, x_2_2__1, x_2_3__1, x_out_2__2, x_out_2__3, z_2_1, z_2_2])

julia> factor(gb[16])
1 * (z_1_2 + 7) * z_1_2
julia> hash(arion_hash.field(-6), arion_hash) # check solution for x_in__1_1
2
julia> hash(arion_hash.field(-2), arion_hash) # check solution for x_in__2_1
2
```
**Usage of Arion Gröbner basis computation experiment** 
```shell
$ cd /Path_to_Arion_implementation/
$ JULIA_NUM_THREADS=N julia ./Arion_groebner_basis_computation.jl # N...number of threads you want to commit to Gröbner basis computation, the higher the faster.
```

**Usage of ArionHash Gröbner basis computation experiment** 
```shell
$ cd /Path_to_Arion_implementation/
$ JULIA_NUM_THREADS=N julia ./ArionHash_groebner_basis_computation.jl # N...number of threads you want to commit to Gröbner basis computation, the higher the faster.
```

**Usage of Arion quotient space dimension evolution experiment**
```shell
$ cd /Path_to_Arion_implementation/
$ JULIA_NUM_THREADS=N julia ./Arion_quotient_space_dimension_evolution.jl # N...number of threads you want to commit to Gröbner basis computation, the higher the faster.
```

**Usage of ArionHash quotient space dimension evolution experiment**
```shell
$ cd /Path_to_Arion_implementation/
$ JULIA_NUM_THREADS=N julia ./ArionHash_quotient_space_dimension_evolution.jl # N...number of threads you want to commit to Gröbner basis computation, the higher the faster.
```
