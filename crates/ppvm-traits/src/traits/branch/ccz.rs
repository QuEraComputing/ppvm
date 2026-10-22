// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// The three-qubit controlled-controlled-Z gate.
pub trait CCZGate<T: Config> {
    /// Apply `CCZ = diag(1, 1, 1, 1, 1, 1, 1, -1)` to distinct qubits.
    ///
    /// The gate is symmetric in its operands and is its own inverse.
    /// If any operand is lost or leaked, the gate has no effect.
    fn ccz(&mut self, a: usize, b: usize, c: usize);

    /// Apply CCZ to each consecutive triple of targets.
    fn ccz_many(&mut self, targets: &[(usize, usize, usize)]) {
        for &(a, b, c) in targets {
            self.ccz(a, b, c);
        }
    }
}
