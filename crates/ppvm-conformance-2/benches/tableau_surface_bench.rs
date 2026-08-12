// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Exhaustive same-binary OLD/NEW microbenchmarks for the public tableau
//! surface. Setup and correctness checks are outside Criterion's timed bodies.

#[path = "tableau_surface/mod.rs"]
mod tableau_surface;

use criterion::{criterion_group, criterion_main};

criterion_group! {
    name = benches;
    config = tableau_surface::criterion_config();
    targets =
        tableau_surface::construction::bench,
        tableau_surface::display::bench,
        tableau_surface::clifford::bench,
        tableau_surface::rotation::bench,
        tableau_surface::measurement::bench,
        tableau_surface::noise::bench,
        tableau_surface::observation::bench,
        tableau_surface::projection::bench,
        tableau_surface::sparse_amplitudes::bench,
}
criterion_main!(benches);
