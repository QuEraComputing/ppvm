// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Complete, operation-by-operation old/new packed-word surface benchmark.

#[path = "word_surface/mod.rs"]
mod word_surface;

use criterion::{criterion_group, criterion_main};

criterion_group!(benches, word_surface::bench);
criterion_main!(benches);
