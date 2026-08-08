// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_conformance_2::{seeded_rng, sym::*};
use ppvm_sym::{Prod as OldProd, Sum as OldSum, Term as OldTerm};
use ppvm_sym_2::{Prod as NewProd, Sum as NewSum, Term as NewTerm};

use super::{assert_real, new_sum, old_sum};

pub(super) fn bench(c: &mut Criterion) {
    bench_eval(c);
    bench_trace(c);
}

fn bench_eval(c: &mut Criterion) {
    let vals = [0.3; 16];
    let mut old_prod = OldProd::sin(0);
    let mut new_prod = NewProd::sin(0);
    for var in 1..16 {
        old_prod.mul_sin(var);
        old_prod.mul_cos(var);
        new_prod.mul_sin(var);
        new_prod.mul_cos(var);
    }

    let mut old_sum = OldSum::new();
    let mut new_sum = NewSum::new();
    for var in 0..16 {
        old_sum.add_term(OldProd::sin(var), 1.0, usize::MAX, f64::EPSILON);
        new_sum.add_term(NewProd::sin(var), 1.0, usize::MAX, f64::EPSILON);
    }
    let old_term = (0..16).fold(OldTerm::from(1.0), |t, var| {
        t + OldTerm::var(var).sin() + OldTerm::var(var).cos()
    });
    let new_term = (0..16).fold(NewTerm::from(1.0), |t, var| {
        t + NewTerm::var(var).sin() + NewTerm::var(var).cos()
    });

    assert_real(old_prod.eval(&vals).unwrap(), new_prod.eval(&vals).unwrap());
    assert_real(old_sum.eval(&vals).unwrap(), new_sum.eval(&vals).unwrap());
    assert_real(old_term.eval(&vals).unwrap(), new_term.eval(&vals).unwrap());

    let mut group = c.benchmark_group("sym/surface/eval");
    group.bench_function("new/prod", |b| b.iter(|| new_prod.eval(&vals).unwrap()));
    group.bench_function("old/prod", |b| b.iter(|| old_prod.eval(&vals).unwrap()));
    group.bench_function("new/sum", |b| b.iter(|| new_sum.eval(&vals).unwrap()));
    group.bench_function("old/sum", |b| b.iter(|| old_sum.eval(&vals).unwrap()));
    group.bench_function("new/term", |b| b.iter(|| new_term.eval(&vals).unwrap()));
    group.bench_function("old/term", |b| b.iter(|| old_term.eval(&vals).unwrap()));
    group.finish();
}

fn bench_trace(c: &mut Criterion) {
    let n = 8;
    let circuit = random_sym_circuit(&mut seeded_rng(0x5A_FACE), n, 120, 8);
    let mut old = old_sum(n);
    old += ("ZIIIIIII", old_seed_coeff(3, 1e-12));
    replay_old(&mut old, &circuit);
    let mut new = new_sum(n);
    new += (NewSymKey::from("ZIIIIIII"), new_seed_coeff(3, 1e-12));
    replay_new(&mut new, &circuit);
    assert_eq!(old.capacity(), new.capacity());

    let old_pattern = ppvm_pauli_word::pattern::PauliPattern::from("Z?*");
    let new_pattern = ppvm_pauli_sum_2::PauliPattern::zero_state();
    let old_trace = ppvm_traits::traits::Trace::trace(&old, &old_pattern);
    let new_trace = ppvm_traits_2::Trace::trace(&new, &new_pattern);
    assert_eq!(
        old_view(&old_trace).monomials,
        new_view(&new_trace).monomials
    );
    assert_real(
        old_trace.eval(&[0.3; 8]).unwrap(),
        new_trace.eval(&[0.3; 8]).unwrap(),
    );

    let mut group = c.benchmark_group("sym/surface/readout");
    group.bench_function("new/trace", |b| {
        b.iter(|| ppvm_traits_2::Trace::trace(&new, &new_pattern))
    });
    group.bench_function("old/trace", |b| {
        b.iter(|| ppvm_traits::traits::Trace::trace(&old, &old_pattern))
    });
    group.finish();
}
