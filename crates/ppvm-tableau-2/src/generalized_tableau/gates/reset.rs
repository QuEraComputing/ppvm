// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Clifford, Measure, Reset, ResetState};

use super::super::GeneralizedTableau;
use crate::VirtualState;

/// Measure, then flip the qubit if it was `|1⟩`.
impl<V: VirtualState> ResetState for GeneralizedTableau<V> {
    fn reset_state(&mut self, q: usize) {
        if self.measure(q) == Some(true) {
            self.x(q);
        }
    }
}

/// No loss or leakage yet, so a full reset is a state reset.
impl<V: VirtualState> Reset for GeneralizedTableau<V> {
    fn reset(&mut self, q: usize) {
        self.reset_state(q);
    }
}
