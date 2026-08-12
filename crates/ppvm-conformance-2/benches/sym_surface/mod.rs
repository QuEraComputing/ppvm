// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Individually named symbolic API benchmarks.
//!
//! Every paired benchmark follows the same discipline:
//! fixtures are built outside the timed body (or in `iter_batched`'s setup
//! closure), and a matching operation is executed once before registration so
//! its output is asserted rather than merely passed to `black_box`.

mod arithmetic;
mod construction;
mod new_only;
mod observability;
mod propagation;
mod readout;

use ppvm_conformance_2::sym::{NewSymSum, OldSymSum};
use ppvm_pauli_sum_2::NoPolicy;

/// Old's exact `NoStrategy` hint. Passing it explicitly to both engines avoids
/// new `NoPolicy`'s intentional `2^14` safety clamp changing benchmark setup.
pub(super) fn matched_capacity(n: usize) -> usize {
    assert!(n > 0 && 2 * n - 1 < usize::BITS as usize);
    1usize << (2 * n - 1)
}

pub(super) fn old_sum(n: usize) -> OldSymSum {
    ppvm_pauli_sum::sum::PauliSum::builder()
        .n_qubits(n)
        .capacity(matched_capacity(n))
        .build()
}

pub(super) fn new_sum(n: usize) -> NewSymSum {
    NewSymSum::with_capacity(n, NoPolicy, matched_capacity(n))
}

pub(super) fn assert_real(old: f64, new: f64) {
    assert!(
        (old - new).abs() < 1e-12,
        "old/new symbolic outputs differ: {old} vs {new}"
    );
}

pub fn bench(c: &mut criterion::Criterion) {
    construction::bench(c);
    arithmetic::bench(c);
    observability::bench(c);
    readout::bench(c);
    propagation::bench(c);
    new_only::bench(c);
}
