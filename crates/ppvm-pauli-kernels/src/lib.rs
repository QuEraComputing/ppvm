// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Bit-level Pauli kernels shared by every PPVM container that conjugates or
//! multiplies Pauli operators.
//!
//! - [`heisenberg`]: one Boolean update rule per Clifford gate, in the
//!   Heisenberg picture `P ↦ G† P G` (so `heisenberg::s` computes `S† P S`).
//!   A rule maps the pre-gate X/Z bits of a Pauli to its post-gate bits and a
//!   sign-flip mask. It is generic over [`Bits`], so the same formula updates
//!   one Pauli (`bool`) or 64 at once (`u64`, one per bit position).
//! - [`row_product`] / [`row_product_phase`]: the packed Pauli-word product and
//!   the exponent `k` of the phase `iᵏ` it emits, built on the per-word
//!   [`word_product`].
//!
//! The crate has no PPVM dependencies; containers wire the rules into their own
//! loops.
//!
//! # Verification
//!
//! The formulas are proved in `lean/PPVM/Kernels/` against exact Gaussian-integer
//! matrices: each rule equals `G† P G` for its gate, and the product phase
//! equals the matrix product's. The tests in `tests/` compare the Rust code
//! against tables generated from that matrix model. `lean/README.md` lists
//! what is not covered by the proofs.

mod bits;
pub mod heisenberg;
mod product;

pub use bits::Bits;
pub use product::{row_product, row_product_phase, word_product};
