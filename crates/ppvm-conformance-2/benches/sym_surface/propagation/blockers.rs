// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Compile-time API blockers, named as such in Criterion output.
//!
//! `Projection` is implemented for the new engine only when its coefficient is
//! `Halvable + One`. `ppvm_sym_2::Term` deliberately has no `Halvable` impl
//! (old's `Term::half()` violates the trait law), so `NewSymSum::p0/p1` do not
//! compile and cannot have an honest paired benchmark. The old operation remains
//! measured under `sym/api_blocker`, never under the comparable surface.
//!
//! `PauliErrorAll` is another true API absence: neither symbolic sum engine
//! implements it. There is therefore no callable operation to register.
//!
//! `AmplitudeDamping` exists on both engines only for `Coeff: num::Float`.
//! Neither old nor new symbolic `Term` implements `Float`, so amplitude damping
//! is likewise unavailable on symbolic sums despite the generic engine impls.

use criterion::{BatchSize, Criterion};
use ppvm_traits::traits::Projection;

use super::fixture;

pub(super) fn bench(c: &mut Criterion) {
    let (old, _) = fixture();
    let mut p0 = old.clone();
    Projection::p0(&mut p0, 0);
    assert!(!p0.data().is_empty());
    let mut p1 = old.clone();
    Projection::p1(&mut p1, 0);
    assert!(!p1.data().is_empty());

    let mut group = c.benchmark_group("sym/api_blocker/projection_new_term_not_halvable");
    group.bench_function("old_only/p0", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                Projection::p0(&mut sum, 0);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.bench_function("old_only/p1", |b| {
        b.iter_batched(
            || old.clone(),
            |mut sum| {
                Projection::p1(&mut sum, 0);
                sum
            },
            BatchSize::SmallInput,
        )
    });
    group.finish();
}
