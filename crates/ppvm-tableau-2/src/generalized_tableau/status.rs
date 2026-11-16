// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Which qubits are lost or leaked. Both are inactive: outside the
//! computational subspace, so gates skip them.

use crate::QubitMask;

/// Inactive and lost qubits as masks; leaked = inactive and not lost.
#[derive(Clone, Debug)]
pub(super) struct Status {
    inactive: QubitMask,
    lost: QubitMask,
    /// Number of inactive qubits, checked first so a loss-free circuit pays
    /// one predictable branch per gate.
    count: usize,
}

impl Status {
    pub(super) fn new(n: usize) -> Self {
        Self {
            inactive: QubitMask::new(n),
            lost: QubitMask::new(n),
            count: 0,
        }
    }

    /// Whether any qubit is inactive.
    pub(super) fn any(&self) -> bool {
        self.count != 0
    }

    pub(super) fn is_inactive(&self, q: usize) -> bool {
        self.any() && self.inactive.contains(q)
    }

    pub(super) fn is_lost(&self, q: usize) -> bool {
        self.any() && self.lost.contains(q)
    }

    pub(super) fn is_leaked(&self, q: usize) -> bool {
        self.is_inactive(q) && !self.lost.contains(q)
    }

    pub(super) fn set_live(&mut self, q: usize) {
        self.mark(q, false, false);
    }

    /// Mark `q` lost; overwrites leakage.
    pub(super) fn set_lost(&mut self, q: usize) {
        self.mark(q, true, true);
    }

    fn mark(&mut self, q: usize, inactive: bool, lost: bool) {
        self.count = self.count - usize::from(self.inactive.contains(q)) + usize::from(inactive);
        self.inactive.set(q, inactive);
        self.lost.set(q, lost);
    }
}
