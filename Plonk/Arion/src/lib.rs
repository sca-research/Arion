#![no_std]

#[cfg(feature = "alloc")]
extern crate alloc;

mod constants;

/// Maximum input width of Arion

/// Arion input width and rounds
pub const MAX_WIDTH: usize = 8;
pub const MAX_ROUNDS: usize = 6;
//pub const WIDTH: usize = 3;
//pub const ROUNDS: usize = 6;
//pub const WIDTH: usize = 4;
//pub const ROUNDS: usize = 5;
pub const WIDTH: usize = 5;
pub const ROUNDS: usize = 5;
//pub const WIDTH: usize = 6;
//pub const ROUNDS: usize = 5;
//pub const WIDTH: usize = 7;
//pub const ROUNDS: usize = 5;
//pub const WIDTH: usize = 8;
//pub const ROUNDS: usize = 4;

/// Strategies implemented for the Arion algorithm.
pub mod strategies;

/// Module containing a fixed-length Arion hash implementation
pub mod perm_uses;

/// Reference implementation for the Arion Sponge hash function
pub mod sponge;

/// The module handling poseidon-trees.
#[cfg(feature = "alloc")]
pub mod tree;
