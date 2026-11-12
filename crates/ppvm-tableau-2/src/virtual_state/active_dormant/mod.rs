// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`ActiveDormant`]: `|χ⟩ = |φ⟩ ⊗ |0…0⟩`, with `φ` a vector over the active
//! virtual qubits behind the [`ActiveState`] interface, and every other
//! (dormant) virtual qubit in `|0⟩`.
//!
//! - `axis`: Paulis on the axes ([`AxisPauli`]) and the collapse rotation
//!   ([`Rotation`]) shared with the tableau.
//! - `state`: the [`ActiveState`] trait.
//! - `dense`: [`DenseState`], an exact vector of `2^k` amplitudes.

mod axis;
mod dense;
mod state;

use num_complex::Complex64;

use super::VirtualState;
use crate::VirtualView;

pub use axis::{AxisPauli, Rotation};
pub use dense::DenseState;
pub use state::ActiveState;

/// `|χ⟩ = |φ⟩ ⊗ |0…0⟩`: amplitudes only for the active virtual qubits.
///
/// Non-Clifford gates activate virtual qubits (localization), and measurements
/// retire them; with no active qubits it is the stabilizer state `C |0…0⟩`.
#[derive(Clone, Debug, PartialEq)]
pub struct ActiveDormant<S: ActiveState> {
    /// Virtual qubit of each axis of `φ`: axis `j` is qubit `active[j]`.
    active: Vec<usize>,
    phi: S,
}

impl<S: ActiveState> ActiveDormant<S> {
    /// The virtual qubit of each axis of `φ`.
    pub fn active(&self) -> &[usize] {
        &self.active
    }

    pub fn phi(&self) -> &S {
        &self.phi
    }
}

impl<S: ActiveState + Default> VirtualState for ActiveDormant<S> {
    fn new(_num_qubits: usize) -> Self {
        Self {
            active: Vec::new(),
            phi: S::default(),
        }
    }

    fn rotate(
        &mut self,
        _t: &mut VirtualView<'_>,
        _qubits: &[usize],
        _a: Complex64,
        _b: Complex64,
    ) {
        todo!("localization and the Pauli rotation on φ (step 7c)")
    }

    /// Deterministic if `O` has no X/Y, else 50/50 (a dormant pivot).
    fn z_weights(&mut self, t: &mut VirtualView<'_>, q: usize) -> [f64; 2] {
        debug_assert_eq!(self.active.len(), self.phi.num_axes());
        if !self.active.is_empty() {
            todo!("measurement with active qubits (step 7d)");
        }
        let o = t.iz(q);
        if o.x.iter().all(|&w| w == 0) {
            if o.negative { [0.0, 1.0] } else { [1.0, 0.0] }
        } else {
            [0.5, 0.5]
        }
    }

    fn collapse(&mut self, t: &mut VirtualView<'_>, q: usize, m: bool) {
        if !self.active.is_empty() {
            todo!("measurement with active qubits (step 7d)");
        }
        let outcome = t.measure(q, || m);
        debug_assert_eq!(outcome, m, "collapsing onto an impossible outcome");
    }

    fn truncation_error(&self) -> f64 {
        self.phi.truncation_error()
    }
}
