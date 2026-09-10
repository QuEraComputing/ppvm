// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Measurement, projection, and reset operations.

use super::{Clifford, CliffordExtensions};

/// Loss-aware projective computational-basis measurement.
///
pub trait Measure {
    /// Measure `qubit`; `None` if the qubit has been lost.
    fn measure<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R) -> Option<bool>;

    /// Measure each target in order, one result per target.
    fn measure_many<R: rand::Rng + ?Sized>(
        &mut self,
        targets: &[usize],
        rng: &mut R,
    ) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q, rng)).collect()
    }
}

// Reset one qubit to a computational/Pauli basis state.
pub trait Reset: Clifford + CliffordExtensions {
    /// Reset one qubit to `|0⟩` (stim `R`/`RZ`).
    fn reset<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R);

    /// stim `RZ` alias — reset to `|0⟩`.
    fn reset_z<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R) {
        self.reset(qubit, rng)
    }

    /// stim `RX` — reset to `|+⟩`.
    fn reset_x<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R) {
        self.reset(qubit, rng);
        self.h(qubit);
    }

    /// stim `RY` — reset to `|i⟩`.
    fn reset_y<R: rand::Rng + ?Sized>(&mut self, qubit: usize, rng: &mut R) {
        self.reset(qubit, rng);
        self.h(qubit);
        self.s(qubit);
    }

    /// Explicit batched reset to `|0⟩`.
    fn reset_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], rng: &mut R) {
        for &q in targets {
            self.reset(q, rng);
        }
    }

    /// Explicit batched `RZ` alias.
    fn reset_z_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], rng: &mut R) {
        self.reset_many(targets, rng)
    }

    /// Explicit batched `RX`.
    fn reset_x_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], rng: &mut R) {
        for &q in targets {
            self.reset_x(q, rng);
        }
    }

    /// Explicit batched `RY`.
    fn reset_y_many<R: rand::Rng + ?Sized>(&mut self, targets: &[usize], rng: &mut R) {
        for &q in targets {
            self.reset_y(q, rng);
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
