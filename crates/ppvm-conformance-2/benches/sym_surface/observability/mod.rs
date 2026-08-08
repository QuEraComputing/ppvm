// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clone and observable trait surface for each symbolic representation.

mod blockers;
mod product;
mod sum;
mod term;

pub(super) fn bench(c: &mut criterion::Criterion) {
    product::bench(c);
    sum::bench(c);
    term::bench(c);
}
