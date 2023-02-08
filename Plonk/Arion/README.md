# ArionHash Plonk Reference Implementation

Contains a pure [Rust](https://www.rust-lang.org/) reference implementation of ArionHash in [Dusk Network Plonk](https://github.com/dusk-network/plonk).

It is based on the Dusk Network [Hades252](https://github.com/dusk-network/Hades252) and [Poseidon252](https://github.com/dusk-network/Poseidon252).

The code for generating Merkle trees and benchmarking is forked from [Poseidon252](https://github.com/dusk-network/Poseidon252) (commit `fc8e47639bb036967d7e6edf3b4d581eb60491cf`).

## Usage

- The [constant_generation](constant_generation) subdirectory contains a [SageMath](https://www.sagemath.org/) Jupyter Notebook to generate constants for Arion in the `BLS12-381` prime field.
To instantiate new round constants, run all cells of the Notebook and copy the outputs into [constants.rs](src/constants.rs).

- With and rounds of ArionHash are hardcoded, to modify them edit [lib.rs](src/lib.rs).

- To benchmark Arion execute the command
```bash
cargo bench
```
