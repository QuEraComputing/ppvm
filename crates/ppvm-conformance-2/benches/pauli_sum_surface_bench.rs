// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Complete old/new Criterion surface comparison for the Pauli-sum engines.
//! Every paired case uses `[u8; 8]`, equal capacity, equal policy semantics, and
//! identical prepared support. Setup cloning is outside timing except where clone
//! or construction itself is the target.

#[path = "pauli_sum_surface/mod.rs"]
mod pauli_sum_surface;

use criterion::{criterion_group, criterion_main};

criterion_group!(
    benches,
    pauli_sum_surface::construction::bench,
    pauli_sum_surface::inspection::bench,
    pauli_sum_surface::clifford::bench,
    pauli_sum_surface::rotation_one::bench,
    pauli_sum_surface::rotation_two::bench,
    pauli_sum_surface::noise::bench,
    pauli_sum_surface::truncation::bench,
    pauli_sum_surface::loss::bench,
    pauli_sum_surface::algebra::bench,
    pauli_sum_surface::projection::bench,
    pauli_sum_surface::representation::bench,
    pauli_sum_surface::new_only::bench,
);
criterion_main!(benches);
