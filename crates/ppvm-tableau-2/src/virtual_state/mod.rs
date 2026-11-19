// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The state `|χ⟩` of the virtual qubits, with `|ψ⟩ = C |χ⟩`, behind the
//! [`VirtualState`] interface.
//!
//! - `active_dormant`: [`ActiveDormant`], `|χ⟩ = |φ⟩ ⊗ |0…0⟩` with `φ` over
//!   the active qubits, and the storage of `φ`.

mod active_dormant;

use num_complex::Complex64;

use crate::VirtualView;

pub use active_dormant::{ActiveDormant, ActiveState, AxisPauli, DenseState, RotationCollapse};

/// The `2ⁿ`-dimensional state `|χ⟩` of the `n` virtual qubits, in a
/// representation of the implementor's choice.
///
/// It sees the tableau only through a [`VirtualView`]: it may change basis
/// (`C ← C·V` with `χ ← V† χ`, so `|ψ⟩` is unchanged), but never applies a
/// physical gate. Arguments are physical qubits; the implementor reads
/// `O = C⁻¹ Z_q C` from the view. The caller samples outcomes, so weights are
/// unnormalized and `rotate` need not be unitary.
pub trait VirtualState: Clone {
    /// `|0…0⟩` on `num_qubits` virtual qubits.
    fn new(num_qubits: usize) -> Self;

    /// Whether `χ = |0…0⟩`, so `|ψ⟩ = C |0…0⟩` is a stabilizer state and the
    /// caller may measure on the tableau alone. `false` is always correct.
    fn is_zero(&self) -> bool {
        false
    }

    /// `χ ← (a·I + b·O) χ` with `O = C⁻¹ (∏ Z_q) C` over `qubits`.
    fn rotate(&mut self, t: &mut VirtualView<'_>, qubits: &[usize], a: Complex64, b: Complex64);

    /// Unnormalized weights of the outcomes `0` and `1` of `Z_q`. Takes `&mut`
    /// so an implementor may change basis here and reuse it in `collapse`.
    fn z_weights(&mut self, t: &mut VirtualView<'_>, q: usize) -> [f64; 2];

    /// Project onto outcome `m` of `Z_q` and renormalize. The outcome must have
    /// nonzero weight.
    fn collapse(&mut self, t: &mut VirtualView<'_>, q: usize, m: bool);
}
