// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[path = "mixture_surface_bench.rs"]
mod mixture_surface_bench;

use criterion::{criterion_group, criterion_main};

criterion_group!(benches, mixture_surface_bench::benchmark);
criterion_main!(benches);
