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

    /// Localization: if `O` has X/Y on dormant qubits, CX from the lowest of
    /// them (in `|0⟩`, so `χ` is unchanged) clears the others, and it becomes
    /// active. Then `O` acts as `±P` on the active axes and as `+1` elsewhere.
    fn rotate(&mut self, t: &mut VirtualView<'_>, qubits: &[usize], a: Complex64, b: Complex64) {
        debug_assert_eq!(self.active.len(), self.phi.num_axes());
        let mut o = t.iz_product(qubits);
        let dormant_x: Vec<usize> = o
            .row()
            .x_support()
            .filter(|v| !self.active.contains(v))
            .collect();
        if let Some((&pivot, rest)) = dormant_x.split_first() {
            t.fan_out(pivot, rest, &[], |x, z| ([x, z], false));
            assert!(self.active.len() < 64, "at most 64 active axes");
            self.active.push(pivot);
            self.phi.activate();
            o = t.iz_product(qubits);
        }
        let o = o.row();
        let mut p = AxisPauli { x: 0, z: 0 };
        for (j, &v) in self.active.iter().enumerate() {
            p.x |= u64::from(o.x_bit(v)) << j;
            p.z |= u64::from(o.z_bit(v)) << j;
        }
        debug_assert_eq!(o.x_support().count(), p.x.count_ones() as usize);
        self.phi.rotate(p, a, if o.negative { -b } else { b });
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
