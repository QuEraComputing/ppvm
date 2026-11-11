// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The active part `φ` of the state `C (|φ⟩ ⊗ |0…0⟩)`: a vector over the
//! active axes `0..k`, behind the [`ActiveState`] interface.
//!
//! - `axis`: Paulis on the axes ([`AxisPauli`]) and the collapse rotation
//!   ([`Rotation`]) shared with the tableau.
//! - `state`: the [`ActiveState`] trait.
//! - `dense`: [`DenseState`], an exact vector of `2^k` amplitudes.

mod axis;
mod dense;
mod state;

pub use axis::{AxisPauli, Rotation};
pub use dense::DenseState;
pub use state::ActiveState;
