// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Reset qubits to a computational/Pauli basis state, mirroring stim's
/// `R`/`RZ`/`RX`/`RY` reset family with broadcasting.
pub trait Reset: crate::traits::Clifford + crate::traits::CliffordExtensions {
    /// Reset each target to `|0⟩` (stim `R`/`RZ`).
    fn reset(&mut self, targets: impl crate::traits::Targets);
    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z(&mut self, targets: impl crate::traits::Targets) {
        self.reset(targets)
    }
    /// stim `RX` — reset to `|+⟩`.
    fn reset_x(&mut self, targets: impl crate::traits::Targets) {
        let qs: Vec<usize> = targets.each().collect();
        self.reset(qs.as_slice());
        self.h(qs.as_slice());
    }
    /// stim `RY` — reset to `|i⟩`.
    fn reset_y(&mut self, targets: impl crate::traits::Targets) {
        let qs: Vec<usize> = targets.each().collect();
        self.reset(qs.as_slice());
        self.h(qs.as_slice());
        self.s(qs.as_slice());
    }
}
