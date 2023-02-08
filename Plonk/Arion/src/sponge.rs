mod hash;

#[cfg(feature = "alloc")]
mod gadget;

pub mod truncated;

pub use hash::hash;

#[cfg(feature = "alloc")]
pub use gadget::gadget;
