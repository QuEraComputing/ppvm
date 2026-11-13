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
//! - `physical_gates`: physical gates `C ← G·C`, which mix rows.
//! - `virtual_gates`: changes of basis `C ← C·V`, which rewrite columns.
//! - `measure`: Z measurement and reset of the stabilizer state.
//! - `transposed`: a column-major copy for batches of changes of basis.
//! - `view`: [`VirtualView`], the virtual side handed to a virtual state.
//!
//! Where rows sit in the arena is decided only in `data`: by `offsets` (used
//! by `x_at`, `z_at`) and by `for_each_row`, which walks the same blocks.

mod bitword;
mod data;
mod measure;
mod physical_gates;
#[cfg(test)]
mod test_support;
mod transposed;
mod view;
mod virtual_gates;

pub use data::{InverseTableau, RowRef};
pub use view::{OwnedRow, VirtualView};
