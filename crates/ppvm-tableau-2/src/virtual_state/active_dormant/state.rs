// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The interface of the active part.

use num_complex::Complex64;

use super::axis::{AxisPauli, Rotation};

/// Storage for `φ` over the active axes `0..k`. Only three things happen to
/// it: a Pauli rotation, an activation, and a collapse; every other change of
/// basis is trivial on the state. It sees axes only, never the tableau.
///
/// `P` is unsigned: the caller folds signs into the coefficients and the
/// outcome. Weights are unnormalized, so `rotate` need not be unitary.
pub trait ActiveState: Clone {
    fn num_axes(&self) -> usize;

    /// Append axis `k` in `|0⟩`.
    fn activate(&mut self);

    /// `φ ← a·φ + b·P φ`.
    fn rotate(&mut self, p: AxisPauli, a: Complex64, b: Complex64);

    /// Unnormalized weights `‖½(1 + P) φ‖²` and `‖½(1 − P) φ‖²` of the
    /// outcomes `+1` and `−1`.
    fn branch_weights(&self, p: AxisPauli) -> [f64; 2];

    /// Project onto outcome `m` of `p` (eigenvalue `(−1)^m`), apply `V†` for
    /// `v = Rotation::collapse(p, m)`, drop axis `v.pivot` (now `|0⟩`), and
    /// renormalize. Axes above the pivot shift down by one. The branch must
    /// have nonzero weight.
    fn collapse(&mut self, p: AxisPauli, m: bool, v: &Rotation);
}
