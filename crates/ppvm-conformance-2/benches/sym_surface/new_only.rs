// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Operations intentionally lacking an old benchmark twin.

use criterion::Criterion;
use ppvm_sym_2::{GaussianInt, Prod, Sum, Term};

pub(super) fn bench(c: &mut Criterion) {
    bench_eval_complex(c);
    bench_exact_ring(c);
}

fn bench_eval_complex(c: &mut Criterion) {
    let vals = [0.3; 8];
    let mut prod = Prod::sin(0);
    for var in 1..8 {
        prod.mul_sin(var);
        prod.mul_cos(var);
    }
    prod.add_phase(1);

    let mut sum = Sum::new();
    for var in 0..8 {
        let mut p = Prod::sin(var);
        p.add_phase((var % 4) as u8);
        sum.add_term(p, 1.0, usize::MAX, f64::EPSILON);
    }
    let term = (0..8).fold(Term::from(1.0), |t, var| {
        t + Term::var(var).sin().mul_phase((var % 4) as u8)
    });

    let prod_value = prod.eval_complex(&vals).unwrap();
    let sum_value = sum.eval_complex(&vals).unwrap();
    let term_value = term.eval_complex(&vals).unwrap();
    assert!(prod_value.re.is_finite() && prod_value.im.is_finite());
    assert!(sum_value.re.is_finite() && sum_value.im.is_finite());
    assert!(term_value.re.is_finite() && term_value.im.is_finite());

    let mut group = c.benchmark_group("sym/new_only/eval_complex");
    group.bench_function("new/prod", |b| b.iter(|| prod.eval_complex(&vals).unwrap()));
    group.bench_function("new/sum", |b| b.iter(|| sum.eval_complex(&vals).unwrap()));
    group.bench_function("new/term", |b| b.iter(|| term.eval_complex(&vals).unwrap()));
    group.finish();
}

fn bench_exact_ring(c: &mut Criterion) {
    let lhs = GaussianInt::new(17, -9);
    let rhs = GaussianInt::new(-4, 11);
    assert_eq!(GaussianInt::from_int(17), GaussianInt::new(17, 0));
    assert_eq!(lhs + rhs, GaussianInt::new(13, 2));
    assert_eq!(lhs * rhs, GaussianInt::new(31, 223));
    assert_eq!(lhs.norm_sq(), 370);

    let mut group = c.benchmark_group("sym/new_only/exact_ring");
    group.bench_function("new/construct", |b| b.iter(|| GaussianInt::new(17, -9)));
    group.bench_function("new/from_int", |b| b.iter(|| GaussianInt::from_int(17)));
    group.bench_function("new/add", |b| b.iter(|| lhs + rhs));
    group.bench_function("new/multiply", |b| b.iter(|| lhs * rhs));
    group.bench_function("new/norm_sq", |b| b.iter(|| lhs.norm_sq()));
    group.finish();
}
