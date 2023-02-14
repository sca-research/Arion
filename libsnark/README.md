# zkp_hash
This library contains implementations, test files and benchmarks of template-parametrized ZK-native permutations and Merkle Trees for ZK verifiable Merkle path authentication.

## Repository structure
The repository is organized as follows:
- `include`: contains headers and templates for the project.
    - `gadget`: contains the code to generate R1CS, and related utilities.
    - `hash`: contains plain implementations of permutation functions.
    - `tree`: contains plain implementations of tree data structures.
    - `util`: contains utility functions.
- `src`: contains benchmark source files and other code which uses the library
- `test`: contains test source files
- `Makefile`: compiles tests and benchmarks.

During compilation, the following additional directories/files are created:
- `bin`: contains executables.
- `build`: contains object files.
- `lib`: contains the static library of the project.

During benchmark execution, the following additional directories/files are created:
- `log`: contains the benchmark logs, in chronological order.

## Requirements
A C++20 compiler and `git`.
To install them run
```shell
sudo apt-get install git-all
```
and
```shell
sudo apt install gcc-11 g++-11
sudo ln -fs /usr/bin/g++-11 /usr/bin/g++
sudo ln -fs /usr/bin/g++-11 /usr/bin/c++
sudo ln -fs /usr/bin/gcc-11 /usr/bin/gcc
```

## Configuring, building and testing
This repository depends on the the following fork+branch of the official `libsnark` library:
> https://github.com/alexander-zw/libsnark/tree/update-libff

If you didn't install `libsnark` manually, then simply run
```shell
sudo ./setup.sh
```
to install `libsnark` and build the library. 

It is possible to configure the project parameters via the `Makefile`, which is self-documented.

To build the project, simply run:
```shell
make
```

To launch all tests, simply run:
```shell
./bin/test.sh
```

You can also execute specific tests by directly running the corresponding binaries in the `bin` directory.

## Benchmarking
For now, to modify the benchmark parameters you have to directly edit the `main()` function body of 
the benchmark files in the `src` directory, and rebuild the project.
Instructions on how to change the code are given as comments in the source files.
