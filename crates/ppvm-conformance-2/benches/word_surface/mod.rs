// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Matched-storage word benchmarks. Setup/reset and correctness checks stay
//! outside timed closures; every operation consumes `black_box`ed inputs.

mod common;
mod construction;
mod lossy;
mod lossy_branch;
mod lossy_gates;
mod lossy_mutation;
mod observation;
mod ordinary;
mod ordinary_gates;
mod ordinary_mutation;
mod pattern;
mod phased;
mod phased_access;
mod phased_gates;

pub fn bench(c: &mut criterion::Criterion) {
    construction::bench(c);
    observation::bench(c);
    ordinary::bench(c);
    ordinary_mutation::bench(c);
    lossy::bench(c);
    lossy_branch::bench(c);
    lossy_mutation::bench(c);
    phased::bench(c);
    phased_access::bench(c);
    ordinary_gates::bench(c);
    lossy_gates::bench(c);
    phased_gates::bench(c);
    pattern::bench(c);
}
