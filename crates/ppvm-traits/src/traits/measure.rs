// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Projective Z-basis measurement returning a bare boolean outcome.
pub trait Measure {
    /// Measure qubit `addr0` in the computational basis. Returns
    /// `true` for outcome `|1⟩`, `false` for `|0⟩`.
    fn measure(&mut self, addr0: usize) -> bool;
}

/// Loss-aware Z-basis measurement.
pub trait LossyMeasure {
    /// Measure qubit `addr0`. Returns `Some(bit)` for an in-subspace
    /// outcome, or `None` if the qubit has been lost.
    fn measure(&mut self, addr0: usize) -> Option<bool>;

    /// Explicit batched measurement: measure each target in order,
    /// returning one result per target.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}
