// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod operator;
mod product;
mod setters;
mod sum;
mod term;

pub(super) fn bench(c: &mut criterion::Criterion) {
    operator::bench(c);
    product::bench(c);
    setters::bench(c);
    sum::bench(c);
    term::bench(c);
}
