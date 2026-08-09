// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::time::Duration;

use criterion::{BatchSize, Criterion};
use ppvm_conformance_2::tableau::Driver;
use ppvm_pauli_sum::config::fx64hash::Byte8F64;

pub mod clifford;
pub mod construction;
pub mod display;
pub mod measurement;
mod measurement_reset;
mod measurement_scratch;
pub mod noise;
pub mod observation;
pub mod projection;
pub mod rotation;
pub mod sparse_amplitudes;
mod sparse_amplitudes_access;

pub type OldBare = ppvm_tableau::data::Tableau<Byte8F64<2>>;
pub type NewBare = ppvm_tableau_2::Tableau<[usize; 2]>;
pub type OldGen = ppvm_conformance_2::tableau::OldWide;
pub type NewGen = ppvm_conformance_2::tableau::NewWide;

pub const WIDTHS: [usize; 6] = [8, 32, 63, 64, 65, 96];
pub const SEED: u64 = 0x5eed_5eed_cafe_babe;
pub const THRESHOLD: f64 = 1e-12;
pub const COEFFICIENT_CAPACITY: usize = 1024;

pub fn criterion_config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(200))
        .measurement_time(Duration::from_millis(700))
        .sample_size(20)
}

pub fn bare_pair(n: usize) -> (OldBare, NewBare) {
    (OldBare::new_with_seed(n, SEED), NewBare::new(n))
}

pub fn gen_pair(n: usize) -> (OldGen, NewGen) {
    let mut old: OldGen = Driver::new_seeded(n, THRESHOLD, SEED);
    let mut new: NewGen = Driver::new_seeded(n, THRESHOLD, SEED);
    let old_len = old.coefficients.len();
    ppvm_tableau::sparsevec::SparseVector::reserve(
        &mut old.coefficients,
        COEFFICIENT_CAPACITY - old_len,
    );
    let new_len = new.coefficients.len();
    new.coefficients.reserve(COEFFICIENT_CAPACITY - new_len);
    (old, new)
}

pub fn prepared_bare(n: usize) -> (OldBare, NewBare) {
    let (mut old, mut new) = bare_pair(n);
    for q in (0..n).step_by(3) {
        ppvm_traits::traits::Clifford::h(&mut old, q);
        ppvm_traits_2::Clifford::h(&mut new, q);
    }
    for q in 0..n - 1 {
        ppvm_traits::traits::Clifford::cnot(&mut old, q, q + 1);
        ppvm_traits_2::Clifford::cnot(&mut new, q, q + 1);
    }
    assert_bare_eq(&old, &new);
    (old, new)
}

pub fn prepared_gen(n: usize) -> (OldGen, NewGen) {
    let (mut old, mut new) = gen_pair(n);
    for q in 0..n.min(8) {
        old.h(q);
        new.h(q);
        old.t(q);
        new.t(q);
    }
    for q in 0..n - 1 {
        old.cnot(q, q + 1);
        new.cnot(q, q + 1);
    }
    assert_gen_eq(&old, &new);
    (old, new)
}

pub fn assert_bare_eq(old: &OldBare, new: &NewBare) {
    assert_eq!(old.to_string(), new.to_string());
}

pub fn assert_gen_eq(old: &OldGen, new: &NewGen) {
    assert_eq!(old.rows(), new.rows());
    assert_eq!(old.record(), new.record());
    assert_eq!(old.lost(), new.lost());
    let (a, b) = (old.coeffs_sorted(), new.coeffs_sorted());
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(x.0, y.0);
        assert!((x.1 - y.1).norm() <= 1e-10);
    }
}

macro_rules! bench_mut_pair {
    ($group:expr, $name:expr, $old_base:expr, $new_base:expr, $old_op:expr, $new_op:expr) => {{
        $group.bench_function(format!("{}/old", $name), |b| {
            b.iter_batched_ref(
                || ($old_base).clone(),
                |value| std::hint::black_box(($old_op)(value)),
                BatchSize::SmallInput,
            )
        });
        $group.bench_function(format!("{}/new", $name), |b| {
            b.iter_batched_ref(
                || ($new_base).clone(),
                |value| std::hint::black_box(($new_op)(value)),
                BatchSize::SmallInput,
            )
        });
    }};
}

pub(crate) use bench_mut_pair;
