# Arion Gröbner Basis Analysis
Contains the [OSCAR](https://oscar.computeralgebra.de/) Gröbner basis attack implementation of Arion and ArionHash.

## Installation
**Requirements**: All scripts have been developed and tested with [Julia](https://julialang.org/) 1.8.4 and [OSCAR](https://oscar.computeralgebra.de/) 0.11.2.

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
julia> Pkg.add(name="Oscar", version="0.11.2")
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
[1]
[5]
[3]
julia> polys = generate_Arion_polynomials(arion=arion,
                                          plain=plain,
                                          cipher=cipher,
                                          field_equations=true);
julia> gb = f4(ideal(polys), nr_thrds=16, info_level=2)
8-element Vector{gfp_mpoly}:
 y_3 + 5
 y_2 + 10
 y_1 + 5
 z_2 + 15
 x_3__1 + 2
 x_2__1 + 11
 x_1__1 + 15
 z_1 + 3
 julia> polys_naive = generate_Arion_polynomials(arion=arion,
                                                 plain=plain,
                                                 cipher=cipher,
                                                 field_equations=true,
                                                 naive_model=true);
julia> gb_naive = f4(ideal(polys_naive), nr_thrds=16, info_level=2)
8-element Vector{gfp_mpoly}:
 y_3 + 5
 y_2 + 10
 y_1 + 5
 z_2 + 15
 x_3__1 + 2
 x_2__1 + 11
 x_1__1 + 15
 z_1 + 3
```

## Examples
**Usage of Arion polynomial model in OSCAR**
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
julia> gb = f4(ideal(polys), nr_thrds=16, info_level=2)
8-element Vector{gfp_mpoly}:
 x_out__3 + 5
 x_out__2 + 7
 z_2 + 16
 x_3__1 + 16
 x_2__1
 x_1__1 + 12
 z_1 + 6
 x_in__1 + 12
 julia> polys_naive = generate_ArionHash_polynomials(arion_hash=arion_hash,
                                                     hash_val=hash_val,
                                                     field_equations=true,
                                                     naive_model=true);
julia> gb_naive = f4(ideal(polys_naive), nr_thrds=16, info_level=2)
8-element Vector{gfp_mpoly}:
 x_out__3 + 5
 x_out__2 + 7
 z_2 + 16
 x_3__1 + 16
 x_2__1
 x_1__1 + 12
 z_1 + 6
 x_in__1 + 12
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
