// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use smallvec::SmallVec;
use stim_parser::prelude::Target;

use super::StimTableau;

const TARGETS_INLINE: usize = 16;

#[inline]
pub(super) fn qubit(t: Target) -> usize {
    t.as_qubit()
        .expect("non-control gate targets are validated as qubits by validate")
}

pub(super) fn qubits(targets: &[Target]) -> SmallVec<[usize; TARGETS_INLINE]> {
    targets.iter().map(|&t| qubit(t)).collect()
}

pub(super) fn qubit_pairs(targets: &[Target]) -> SmallVec<[(usize, usize); TARGETS_INLINE / 2]> {
    targets
        .chunks_exact(2)
        .map(|p| (qubit(p[0]), qubit(p[1])))
        .collect()
}

pub(super) fn has_record_control(targets: &[Target]) -> bool {
    targets.iter().any(|t| matches!(t, Target::Rec(_)))
}

#[inline]
pub(super) fn record_bit(record: &[Option<bool>], k: usize) -> bool {
    record
        .len()
        .checked_sub(k)
        .and_then(|i| record.get(i).copied().flatten())
        .unwrap_or(false)
}

pub(super) fn measure_reset_z<T: StimTableau>(tab: &mut T, q: usize, noise: f64) -> Option<bool> {
    let true_outcome = tab.measure(q);
    if true_outcome == Some(true) {
        tab.x(q);
    }
    let recorded = true_outcome.map(|b| tab.flip_with_prob(b, noise));
    tab.overwrite_last_measurement_record(recorded);
    recorded
}
