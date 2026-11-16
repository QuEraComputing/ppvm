// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Clifford, Reset, ResetState};

use super::super::GeneralizedTableau;
use crate::VirtualState;

/// Measure, then flip the qubit if it was `|1⟩`. Lost and leaked qubits are
/// left alone.
impl<V: VirtualState> ResetState for GeneralizedTableau<V> {
    fn reset_state(&mut self, q: usize) {
        if !self.status.is_inactive(q) && self.sample_z(q) {
            self.x(q);
        }
    }
}

/// A full reset: a lost or leaked qubit becomes live, then its state is reset.
impl<V: VirtualState> Reset for GeneralizedTableau<V> {
    fn reset(&mut self, q: usize) {
        self.status.set_live(q);
        self.reset_state(q);
    }
}
