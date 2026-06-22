// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Reset one qubit to a computational/Pauli basis state.
pub trait Reset: crate::traits::Clifford + crate::traits::CliffordExtensions {
    /// Reset one qubit to `|0⟩` (stim `R`/`RZ`).
    fn reset(&mut self, addr0: usize);
    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z(&mut self, addr0: usize) {
        self.reset(addr0)
    }
    /// stim `RX` — reset to `|+⟩`.
    fn reset_x(&mut self, addr0: usize) {
        self.reset(addr0);
        self.h(addr0);
    }
    /// stim `RY` — reset to `|i⟩`.
    fn reset_y(&mut self, addr0: usize) {
        self.reset(addr0);
        self.h(addr0);
        self.s(addr0);
    }

    /// Explicit batched reset to `|0⟩`.
    fn reset_batch(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset(q);
        }
    }

    /// Explicit batched `RZ` alias.
    fn reset_z_batch(&mut self, targets: &[usize]) {
        self.reset_batch(targets)
    }

    /// Explicit batched `RX`.
    fn reset_x_batch(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_x(q);
        }
    }

    /// Explicit batched `RY`.
    fn reset_y_batch(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_y(q);
        }
    }
}
