# libsnark Implementation of ArionHash
This library contains implementations, test files and benchmarks of template-parametrized ZK-native permutations and Merkle Trees for ZK verifiable Merkle path authentication.
This implementation is updated up to paper-preprint status, any future additions will be made in the following fork: https://github.com/byt3bit/zkp_hash

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
First install requirements for `libsnark`
```shell
sudo apt-get install build-essential git-all libboost-all-dev cmake libgmp-dev libssl-dev libprocps-dev pkg-config gnuplot-x11 libsodium-dev
```
On Ubuntu 20 one additionally needs to install a `C++20` compliant compiler, like `gcc-11+`: 
```shell
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get install gcc-11 g++-11
sudo ln -fs /usr/bin/g++-11 /usr/bin/g++
sudo ln -fs /usr/bin/g++-11 /usr/bin/c++
sudo ln -fs /usr/bin/gcc-11 /usr/bin/gcc
```

## Configuring, building and testing
This repository depends on the following fork+branch of the official `libsnark` library:
> https://github.com/alexander-zw/libsnark/tree/update-libff

If you didn't install `libsnark` manually, then simply run
```shell
sh ./setup.sh
```
to install `libsnark` and build the library. 

It is possible to configure the project parameters via the `Makefile`, which is self-documented.

To build the project, simply run:
```shell
make
```

To launch all tests, simply run:
```shell
sh ./bin/test.sh
```

You can also execute specific tests by directly running the corresponding binaries in the `bin` directory.

## Benchmarking
For now, to modify the benchmark parameters you have to directly edit the `main()` function body of 
the benchmark files in the `src` directory, and rebuild the project.
Instructions on how to change the code are given as comments in the source files.

To execute a benchmark run
```shell
./bin/benchmark_mtree
```
Log files will be generated in `libsnark/log`.
Make sure that you always run the benchmark from the `libsnark` directory, else log files end up in digital Nirvana.
