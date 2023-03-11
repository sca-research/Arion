# Arion
This repository contains implementations of Arion in [SageMath](https://www.sagemath.org/), [OSCAR](https://oscar.computeralgebra.de/), [libsnark](https://github.com/scipr-lab/libsnark) and [Plonk](https://github.com/dusk-network/plonk).

The Arion preprint is available at https://arxiv.org/abs/2303.04639.

## Repository Structure
The repository is divided in four main directories:
- `libsnark` : contains `C++` vanilla and SNARK implementations of our and competitor designs by leveraging the Groth16 framework. It allows one to test and compare such constructions for Merkle Tree commitment verification.
- `OSCAR` : Contains the `Julia` plain implementation of our design, together with the Grobner basis attack analysis.
- `Plonk` : Contains the `Rust` circuit implementation of Arion based on the Plonk gates
- `SageMath`: Contains the `SageMath` plain implementation of Arion.
