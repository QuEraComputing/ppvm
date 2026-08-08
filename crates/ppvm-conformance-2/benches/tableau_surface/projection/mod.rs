// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use num::complex::Complex64;
use ppvm_conformance_2::tableau::Driver;

use super::*;

mod overlap;
mod project;

type Entries = Vec<(Complex64, u128)>;

#[derive(Clone, Copy)]
struct Decomposition {
    phase: u8,
    stabilizer: u128,
    destabilizer: u128,
}

fn projection_pair() -> (OldGen, NewGen) {
    let (mut old, mut new) = gen_pair(8);
    // Fixed support 2^4 = 16: enough to exercise map pairing without an
    // exponentially growing benchmark workload.
    for q in 1..=4 {
        old.h(q);
        new.h(q);
        old.t(q);
        new.t(q);
    }
    old.h(0);
    new.h(0);
    assert_gen_eq(&old, &new);
    (old, new)
}

fn decomposition(old: &OldGen, new: &NewGen, qubit: usize) -> Decomposition {
    let old_d = old.compute_decomposition(qubit, ppvm_traits::Pauli::Z);
    let new_d = new.compute_decomposition(qubit, ppvm_traits_2::Pauli::Z);
    assert_eq!(old_d, new_d);
    Decomposition {
        phase: old_d.0,
        stabilizer: old_d.1,
        destabilizer: old_d.2,
    }
}

fn entries(old: &OldGen, new: &NewGen) -> (Entries, Entries) {
    let old_entries = old.coefficients.to_vec();
    let new_entries = new.coefficients.entries().to_vec();
    assert_eq!(old_entries, new_entries);
    (old_entries, new_entries)
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/projection");
    let (old, new) = projection_pair();
    project::bench(&mut group, &old, &new);
    overlap::bench(&mut group, &old, &new);
    group.finish();
}
