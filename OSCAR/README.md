# Arion Gröbner Basis Analysis
Contains the [OSCAR](https://oscar.computeralgebra.de/) Gröbner basis attack implementation of Arion and ArionHash.

## Installation
**Requirements**: All scripts have been developed and tested with [Julia](https://julialang.org/) 1.7.2 and [OSCAR](https://oscar.computeralgebra.de/) 0.10.2

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
julia> gb_naive = gb = f4(ideal(polys_naive), nr_thrds=16, info_level=2)
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