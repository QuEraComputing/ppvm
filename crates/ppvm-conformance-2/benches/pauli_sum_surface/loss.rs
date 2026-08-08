// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_lossy_pauli_word_2::LossyPauliWord as NewWord;
use ppvm_pauli_sum::config;
use ppvm_pauli_sum::strategy::MaxLossWeight as OldMaxLoss;
use ppvm_pauli_sum::sum::PauliSum;
use ppvm_pauli_sum_2::{HashMapStore, MaxLossWeight as NewMaxLoss, NoPolicy, Sum};
use ppvm_pauli_word::loss::LossyPauliWord as OldWord;
use ppvm_traits::traits::{
    CorrelatedLossChannel as OldCorrelated, LossChannel as OldLoss, NoStrategy,
    ResetLossChannel as OldReset,
};
use ppvm_traits_2::{
    CorrelatedLossChannel as NewCorrelated, LossChannel as NewLoss, ResetLossChannel as NewReset,
};

use super::{CAPACITY, N};

type OldKey = OldWord<[u8; 8], fxhash::FxBuildHasher>;
type NewKey = NewWord<[u8; 8], fxhash::FxBuildHasher>;
type OldCfg = config::fxhash::Byte<8, f64, NoStrategy, OldKey>;
type OldSum = PauliSum<OldCfg>;
type NewSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;
type OldLossCfg = config::fxhash::Byte<8, f64, OldMaxLoss, OldKey>;
type OldLossSum = PauliSum<OldLossCfg>;
type NewLossSum = Sum<HashMapStore<NewKey, f64>, NewMaxLoss>;

fn data() -> Vec<(String, f64)> {
    const SITES: [char; 5] = ['I', 'X', 'Y', 'Z', 'L'];
    (0..160)
        .map(|i| {
            let mut x = i * 7919 + 11;
            let word = (0..N)
                .map(|_| {
                    let p = SITES[x % 5];
                    x /= 5;
                    p
                })
                .collect();
            (word, 0.25 + (i % 17) as f64 / 8.0)
        })
        .collect()
}

fn old_sum() -> OldSum {
    let mut sum = OldSum::builder()
        .n_qubits(N)
        .strategy(NoStrategy)
        .capacity(CAPACITY)
        .build();
    for (word, coeff) in data() {
        sum += (word.as_str(), coeff);
    }
    sum
}

fn new_sum() -> NewSum {
    let mut sum = NewSum::with_capacity(N, NoPolicy, CAPACITY);
    for (word, coeff) in data() {
        sum += (NewKey::from(word.as_str()), coeff);
    }
    sum
}

fn support_old<T: ppvm_traits::config::Config>(sum: &PauliSum<T>) -> Vec<(String, u64)>
where
    for<'a> T::Map: ppvm_traits::traits::ACMapIter<'a, Item = (&'a T::PauliWordType, &'a f64)>,
    T::Coeff: Into<f64>,
    T::PauliWordType: ToString,
{
    let mut out: Vec<_> = sum
        .iter()
        .map(|(w, c)| (w.to_string(), c.to_bits()))
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn support_new<P: ppvm_pauli_sum_2::Policy<NewKey, f64>>(
    sum: &Sum<HashMapStore<NewKey, f64>, P>,
) -> Vec<(String, u64)> {
    let mut out: Vec<_> = sum
        .iter()
        .map(|(w, c)| (w.to_string(), c.to_bits()))
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn paired<FO, FN>(c: &mut Criterion, name: &str, old_op: FO, new_op: FN)
where
    FO: Fn(&mut OldSum) + Copy,
    FN: Fn(&mut NewSum) + Copy,
{
    let (old, new) = (old_sum(), new_sum());
    let (mut op, mut np) = (old.clone(), new.clone());
    old_op(&mut op);
    new_op(&mut np);
    assert_eq!(support_old(&op), support_new(&np));
    let mut group = c.benchmark_group(format!("pauli_sum_surface/loss/{name}"));
    group.bench_function("old", |b| {
        b.iter_batched_ref(|| old.clone(), old_op, BatchSize::LargeInput)
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(|| new.clone(), new_op, BatchSize::LargeInput)
    });
    group.finish();
}

pub fn bench(c: &mut Criterion) {
    paired(
        c,
        "channel",
        |s| s.loss_channel(3, 0.17),
        |s| s.loss_channel(3, 0.17),
    );
    paired(
        c,
        "reset",
        |s| s.reset_loss_channel(3),
        |s| s.reset_loss_channel(3),
    );
    paired(
        c,
        "correlated",
        |s| s.correlated_loss_channel(2, 5, [0.07, 0.11, 0.19]),
        |s| s.correlated_loss_channel(2, 5, [0.07, 0.11, 0.19]),
    );
    max_loss(c, usize::MAX, "max_loss_weight_disabled");
    max_loss(c, 2, "max_loss_weight_active");
}

fn max_loss(c: &mut Criterion, bound: usize, name: &str) {
    let mut old: OldLossSum = OldLossSum::builder()
        .n_qubits(N)
        .strategy(OldMaxLoss(bound))
        .capacity(CAPACITY)
        .build();
    let mut new = NewLossSum::with_capacity(N, NewMaxLoss(bound), CAPACITY);
    for (word, coeff) in data() {
        old += (word.as_str(), coeff);
        new += (NewKey::from(word.as_str()), coeff);
    }
    let (mut op, mut np) = (old.clone(), new.clone());
    op.truncate();
    np.truncate();
    assert_eq!(support_old(&op), support_new(&np));
    let mut group = c.benchmark_group(format!("pauli_sum_surface/truncate/{name}"));
    group.bench_function("old", |b| {
        b.iter_batched_ref(|| old.clone(), |s| s.truncate(), BatchSize::LargeInput)
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(|| new.clone(), |s| s.truncate(), BatchSize::LargeInput)
    });
    group.finish();
}
