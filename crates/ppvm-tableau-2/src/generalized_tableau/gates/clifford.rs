// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{ApplyGateRule, CliffordBatch};

use super::super::GeneralizedTableau;
use crate::VirtualState;

/// Physical gates `C ← G·C` go to the tableau; `|χ⟩` is untouched.
impl<V: VirtualState> ApplyGateRule for GeneralizedTableau<V> {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.tableau.apply1(q, rule);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        self.tableau.apply2(a, b, rule);
    }

    fn apply_sign1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.tableau.apply_sign1(q, rule);
    }
}

impl<V: VirtualState> CliffordBatch for GeneralizedTableau<V> {}
