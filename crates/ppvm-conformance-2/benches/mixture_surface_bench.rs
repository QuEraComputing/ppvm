// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Criterion surface comparison for the legacy `GeneralizedTableauSum` and
//! the `-2` `GeneralizedTableauMixture`.
//!
//! Cargo integration intentionally lives in the already-registered
//! `tableau_mixture_bench` target because this conformance crate's manifest is
//! frozen during the comparison.

#[path = "mixture_surface_bench/clifford.rs"]
mod clifford;
#[path = "mixture_surface_bench/integration.rs"]
mod integration;
#[path = "mixture_surface_bench/lifecycle.rs"]
mod lifecycle;
#[path = "mixture_surface_bench/measure.rs"]
mod measure;
#[path = "mixture_surface_bench/new_only.rs"]
mod new_only;
#[path = "mixture_surface_bench/noise.rs"]
mod noise;
#[path = "mixture_surface_bench/rotations.rs"]
mod rotations;
#[path = "mixture_surface_bench/sampling.rs"]
mod sampling;
#[path = "mixture_surface_bench/scaling.rs"]
mod scaling;
#[path = "mixture_surface_bench/support.rs"]
mod support;

use criterion::Criterion;

pub fn benchmark(c: &mut Criterion) {
    lifecycle::register(c);
    clifford::register(c);
    rotations::register(c);
    measure::register(c);
    noise::register(c);
    sampling::register(c);
    scaling::register(c);
    integration::register(c);
    new_only::register(c);
}
