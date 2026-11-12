// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`VirtualView`]: the part of the tableau a virtual state may use.

use super::{InverseTableau, RowRef};

/// The virtual side of an [`InverseTableau`]: row reads and changes of basis
/// `C ← C·V`, but no physical gates. A virtual state gets this view, so it can
/// change basis together with its own state while `|ψ⟩` stays the same.
pub struct VirtualView<'a> {
    t: &'a mut InverseTableau,
}

impl InverseTableau {
    pub fn virtual_view(&mut self) -> VirtualView<'_> {
        VirtualView { t: self }
    }
}

impl VirtualView<'_> {
    pub fn num_qubits(&self) -> usize {
        self.t.num_qubits()
    }

    /// Row `ix_q = C⁻¹ X_q C`.
    pub fn ix(&self, q: usize) -> RowRef<'_> {
        self.t.ix(q)
    }

    /// Row `iz_q = C⁻¹ Z_q C`.
    pub fn iz(&self, q: usize) -> RowRef<'_> {
        self.t.iz(q)
    }

    /// See [`InverseTableau::change_basis1`].
    pub fn change_basis1(&mut self, v: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.t.change_basis1(v, rule);
    }

    /// See [`InverseTableau::change_basis2`].
    pub fn change_basis2(
        &mut self,
        a: usize,
        b: usize,
        rule: impl Fn([bool; 4]) -> ([bool; 4], bool),
    ) {
        self.t.change_basis2(a, b, rule);
    }

    /// See [`InverseTableau::fan_out`].
    pub fn fan_out(
        &mut self,
        pivot: usize,
        cx_targets: &[usize],
        cz_targets: &[usize],
        last: impl Fn(bool, bool) -> ([bool; 2], bool),
    ) {
        self.t.fan_out(pivot, cx_targets, cz_targets, last);
    }

    /// Measure `Z_q` with every virtual qubit dormant, i.e. on `C |0…0⟩`. See
    /// [`InverseTableau::measure`]; the collapse is a change of basis only.
    pub fn measure(&mut self, q: usize, coin: impl FnOnce() -> bool) -> bool {
        self.t.measure(q, coin)
    }
}
