//! This crate contains implementations of $`\mathbb{F}_2`$ as well as binary extension fields.

#![deny(missing_docs)]
mod f128b;
mod f2;
mod f64b;
mod f8b;
mod small_binary_fields;
pub use f128b::*;
pub use f2::*;
pub use f64b::*;
pub use f8b::*;
pub use small_binary_fields::*;
