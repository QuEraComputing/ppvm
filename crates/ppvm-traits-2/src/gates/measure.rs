// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Measurement, projection, and reset operations.

use super::Clifford;

/// Loss-aware projective computational-basis measurement.
///
pub trait Measure {
    /// Measure `qubit`; `None` if the qubit has been lost.
    fn measure(&mut self, qubit: usize) -> Option<bool>;

    /// Measure each target in order, one result per target.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}

// Reset one qubit to a computational/Pauli basis state.
pub trait Reset: Clifford {
    /// Reset one qubit to `|0⟩` (stim `R`/`RZ`).
    fn reset(&mut self, qubit: usize);

    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z(&mut self, qubit: usize) {
        self.reset(qubit)
    }

    /// stim `RX` — reset to `|+⟩`.
    fn reset_x(&mut self, qubit: usize) {
        self.reset(qubit);
        self.h(qubit);
    }

    /// stim `RY` — reset to `|i⟩`.
    fn reset_y(&mut self, qubit: usize) {
        self.reset(qubit);
        self.h(qubit);
        self.s(qubit);
    }

    /// Explicit batched reset to `|0⟩`.
    fn reset_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset(q);
        }
    }

    /// Explicit batched `RZ` alias.
    fn reset_z_many(&mut self, targets: &[usize]) {
        self.reset_many(targets)
    }

    /// Explicit batched `RX`.
    fn reset_x_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_x(q);
        }
    }

    /// Explicit batched `RY`.
    fn reset_y_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.reset_y(q);
        }
    }
}

/// Projective Z-basis projectors `|0⟩⟨0|` and `|1⟩⟨1|`
pub trait Projection {
    /// Project `qubit` onto `|0⟩`.
    fn p0(&mut self, qubit: usize);
    /// Project `qubit` onto `|1⟩`.
    fn p1(&mut self, qubit: usize);
}
