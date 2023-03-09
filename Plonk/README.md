# Plonk Implementation of ArionHash

This repository contains the pure [Rust](https://www.rust-lang.org/) implementation of ArionHash in the [Dusk Network Plonk](https://github.com/dusk-network/plonk) library.
It is based on the Dusk Network [Hades252](https://github.com/dusk-network/Hades252) and [Poseidon252](https://github.com/dusk-network/Poseidon252).

The code for generating Merkle trees and benchmarking is forked from [Poseidon252](https://github.com/dusk-network/Poseidon252) (commit `fc8e47639bb036967d7e6edf3b4d581eb60491cf`).

## Usage

To benchmark Arion make sure that [Rust](https://www.rust-lang.org/) is installed on your system.
Then one can compile and run the benchmark using the following commands.

```bash
cd Arion
cargo bench
```
To benchmark [Poseidon252](https://github.com/dusk-network/Poseidon252) use the following commands
```bash
cd "Poseidon Benchmark"
cargo bench
```
