// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod clifford;
mod reset;
mod rot1;
mod rot2;
mod tgate;
mod u3;

macro_rules! impl_generalized_tableau_sum_gate {
    ($name:ident) => {
        fn $name(&mut self, targets: impl ppvm_traits::traits::Targets) {
            let targets: Vec<usize> = targets.each().collect();
            self.entries.for_each_mut(|tab, _p| {
                tab.$name(targets.as_slice());
            });
            // The gate mutates every entry's tableau (or no-ops on a
            // lost qubit, in which case the cached fp is still valid).
            // Conservatively clear all cached fingerprints; they'll be
            // recomputed lazily on the next insert_or_update_batch.
            self.entries.mark_keys_dirty();
        }
    };
}
pub(crate) use impl_generalized_tableau_sum_gate;
