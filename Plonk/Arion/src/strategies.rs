//! This module contains an implementation of the `Arion`
//! strategy algorithm specifically designed to work outside of
//! Rank 1 Constraint Systems (R1CS) or other custom Constraint
//! Systems such as Add/Mul/Custom plonk gate-circuits.
//!
//! The inputs of the permutation function have to be explicitly
//! over the BlsScalar Field of the bls12_381 curve so working over
//! `Fq = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`.

/// Strategy for zero-knowledge plonk circuits
mod gadget;

/// Strategy for scalars
mod scalar;

use dusk_bls12_381::BlsScalar;

pub use gadget::GadgetStrategy;
pub use scalar::ScalarStrategy;

use crate::{ROUNDS, WIDTH, constants::{CONSTANTS_G, CONSTANTS_H, CONSTANTS_AFF}};

/// Defines the Arion algorithm.
pub trait Strategy<T: Clone + Copy> {
    fn quinti_s_box(&mut self, input: &mut T);

    fn gtds(&mut self, input: &mut [T], constants_g: &[[BlsScalar; 2]; WIDTH - 1], constants_h: &[BlsScalar; WIDTH - 1]);

    fn mul_matrix(&mut self, input: &mut [T]);

    fn affine_layer(&mut self, input: &mut [T], constants_aff: &[BlsScalar; WIDTH]);

    fn perm(&mut self, data: &mut [T]) {
        self.mul_matrix(data);
        self.affine_layer(data, &[BlsScalar::zero(); WIDTH]);
        for r in 0..ROUNDS {
            self.gtds(data, &CONSTANTS_G[r], &CONSTANTS_H[r]);
            self.affine_layer(data, &CONSTANTS_AFF[r]);
        }
    }
}
