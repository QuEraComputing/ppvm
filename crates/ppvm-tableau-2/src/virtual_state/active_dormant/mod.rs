// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`ActiveDormant`]: `|χ⟩ = |φ⟩ ⊗ |0…0⟩`, with `φ` a vector over the active
//! virtual qubits behind the [`ActiveState`] interface, and every other
//! (dormant) virtual qubit in `|0⟩`.
//!
//! - `axis`: Paulis on the axes ([`AxisPauli`]) and the collapse rotation
//!   ([`RotationCollapse`]) shared with the tableau.
//! - `state`: the [`ActiveState`] trait.
//! - `dense`: [`DenseState`], an exact vector of `2^k` amplitudes.

mod axis;
mod dense;
mod state;

use num_complex::Complex64;

use ppvm_traits_2::ApplyGateRule;

use super::VirtualState;
use crate::{QubitMask, RowRef, VirtualView};

pub use axis::{AxisPauli, RotationCollapse};
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
    /// The qubits of `active`, for word-wise classification.
    mask: QubitMask,
    phi: S,
}

/// How `O = C⁻¹ Z_q C` acts on `|φ⟩ ⊗ |0…0⟩`.
enum Case {
    /// No X/Y on dormant qubits and none on the active axes: `O = ±1` there.
    Deterministic { outcome: bool },
    /// X/Y on a dormant qubit: 50/50, `φ` untouched.
    DormantPivot,
    /// `O = ±P` with `P` a nontrivial Pauli on the active axes.
    Active { p: AxisPauli, negative: bool },
}

impl<S: ActiveState> ActiveDormant<S> {
    /// The virtual qubit of each axis of `φ`.
    pub fn active(&self) -> &[usize] {
        &self.active
    }

    pub fn phi(&self) -> &S {
        &self.phi
    }

    /// `o` on the active axes; its X bits must all be on active qubits.
    fn axis_pauli(&self, o: RowRef<'_>) -> AxisPauli {
        let mut p = AxisPauli { x: 0, z: 0 };
        for (j, &v) in self.active.iter().enumerate() {
            p.x |= u64::from(o.x_bit(v)) << j;
            p.z |= u64::from(o.z_bit(v)) << j;
        }
        debug_assert_eq!(o.x_support().count(), p.x.count_ones() as usize);
        p
    }

    fn classify(&self, o: RowRef<'_>) -> Case {
        if o.x_outside(&self.mask).next().is_some() {
            return Case::DormantPivot;
        }
        let p = self.axis_pauli(o);
        if p.is_identity() {
            Case::Deterministic {
                outcome: o.negative,
            }
        } else {
            Case::Active {
                p,
                negative: o.negative,
            }
        }
    }
}

/// The view with axis `j` as virtual qubit `active[j]`, so a [`RotationCollapse`]'s
/// gates change the tableau's basis.
struct OnAxes<'v, 'a> {
    t: &'v mut VirtualView<'a>,
    active: &'v [usize],
}

impl ApplyGateRule for OnAxes<'_, '_> {
    type B = bool;

    fn apply1(&mut self, j: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.t.change_basis1(self.active[j], rule);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        self.t.change_basis2(self.active[a], self.active[b], rule);
    }
}

impl<S: ActiveState + Default> VirtualState for ActiveDormant<S> {
    fn new(num_qubits: usize) -> Self {
        Self {
            active: Vec::new(),
            mask: QubitMask::new(num_qubits),
            phi: S::default(),
        }
    }

    fn is_zero(&self) -> bool {
        self.active.is_empty()
    }

    /// Localization: if `O` has X/Y on dormant qubits, CX from the lowest of
    /// them (in `|0⟩`, so `χ` is unchanged) clears the others, and it becomes
    /// active. Then `O` acts as `±P` on the active axes and as `+1` elsewhere.
    fn rotate(&mut self, t: &mut VirtualView<'_>, qubits: &[usize], a: Complex64, b: Complex64) {
        debug_assert_eq!(self.active.len(), self.phi.num_axes());
        let mut o = t.iz_product(qubits);
        let dormant_x: Vec<usize> = o.row().x_outside(&self.mask).collect();
        if let Some((&pivot, rest)) = dormant_x.split_first() {
            t.fan_out(pivot, rest, &[], |x, z| ([x, z], false));
            assert!(self.active.len() < 64, "at most 64 active axes");
            self.active.push(pivot);
            self.mask.set(pivot, true);
            self.phi.activate();
            o = t.iz_product(qubits);
        }
        let o = o.row();
        let p = self.axis_pauli(o);
        self.phi.rotate(p, a, if o.negative { -b } else { b });
    }

    /// Deterministic, 50/50 for a dormant pivot, or the branch weights of `±P`.
    fn z_weights(&mut self, t: &mut VirtualView<'_>, q: usize) -> [f64; 2] {
        debug_assert_eq!(self.active.len(), self.phi.num_axes());
        match self.classify(t.iz(q)) {
            Case::Deterministic { outcome: false } => [1.0, 0.0],
            Case::Deterministic { outcome: true } => [0.0, 1.0],
            Case::DormantPivot => [0.5, 0.5],
            Case::Active { p, negative } => {
                let [plus, minus] = self.phi.branch_weights(p);
                if negative {
                    [minus, plus]
                } else {
                    [plus, minus]
                }
            }
        }
    }

    /// The active case measures `P` with outcome `m ⊕ negative`, changes basis
    /// so that its pivot axis ends in `|0⟩`, and retires that axis.
    fn collapse(&mut self, t: &mut VirtualView<'_>, q: usize, m: bool) {
        match self.classify(t.iz(q)) {
            Case::Deterministic { outcome } => {
                debug_assert_eq!(outcome, m, "collapsing onto an impossible outcome");
            }
            Case::DormantPivot => t.collapse_dormant(q, &self.mask, m),
            Case::Active { p, negative } => {
                let m = m != negative;
                let v = RotationCollapse::collapse(p, m);
                v.apply_gates(&mut OnAxes {
                    t,
                    active: &self.active,
                });
                self.phi.collapse(p, m, &v);
                let retired = self.active.remove(v.pivot);
                self.mask.set(retired, false);
            }
        }
    }
}
