// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The inverse stabilizer tableau `C⁻¹` of a Clifford `C`.
//!
//! Row `ix_q` is `C⁻¹ X_q C` and row `iz_q` is `C⁻¹ Z_q C`, each stored as
//! packed X and Z bits plus one sign bit (rows are Hermitian, so their phase
//! is ±1). Rows are indexed `ix_q = q`, `iz_q = n + q`.
//!
//! - `data`: the storage, row reads, and where rows sit in the arena.
//! - `bitword`: the storage word; the only code that knows its width.
//! - `physical`: physical gates `C ← G·C`, which mix rows.
//! - `virtual_side`: changes of basis `C ← C·V`, which rewrite columns.
//!
//! Where rows sit in the arena is decided only by `offsets` (used by `x_at`,
//! `z_at` and `row_offsets`).

mod bitword;
mod data;
mod physical;
#[cfg(test)]
mod test_support;
mod virtual_side;

pub use data::{InverseTableau, RowRef};
