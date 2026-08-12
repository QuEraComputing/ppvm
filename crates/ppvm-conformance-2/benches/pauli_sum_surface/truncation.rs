// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeff, CombinedStrategy as OldCombined, MaxPauliWeight as OldWeight,
};
use ppvm_pauli_sum::sum::PauliSum;
use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeff, CombinedPolicy as NewCombined, HashMapStore,
    MaxPauliWeight as NewWeight, PauliWord, Sum,
};

use super::{CAPACITY, N, NewKey, OldKey};

type OldCoeffSum = PauliSum<ByteF64<8, OldCoeff>>;
type NewCoeffSum = Sum<HashMapStore<NewKey, f64>, NewCoeff>;
type OldWeightSum = PauliSum<ByteF64<8, OldWeight>>;
type NewWeightSum = Sum<HashMapStore<NewKey, f64>, NewWeight>;
type OldCombinedSum = PauliSum<ByteF64<8, OldCombined<OldCoeff, OldWeight>>>;
type NewCombinedSum = Sum<HashMapStore<NewKey, f64>, NewCombined<NewCoeff, NewWeight>>;

const DATA: [(&str, f64); 8] = [
    ("IIIIIIII", 1e-8),
    ("XIIIIIII", 0.25),
    ("ZZIIIIII", 0.75),
    ("XYZIIIII", 1e-7),
    ("XXXXIIII", 1.25),
    ("ZZZZZIII", 0.5),
    ("XYZXYZII", 2.0),
    ("ZZZZZZZZ", 1e-9),
];

macro_rules! fill_old {
    ($ty:ty, $strategy:expr) => {{
        let mut sum: $ty = PauliSum::builder()
            .n_qubits(N)
            .strategy($strategy)
            .capacity(CAPACITY)
            .build();
        for (word, coeff) in DATA {
            sum += (word, coeff);
        }
        sum
    }};
}

macro_rules! fill_new {
    ($ty:ty, $policy:expr) => {{
        let mut sum: $ty = Sum::with_capacity(N, $policy, CAPACITY);
        for (word, coeff) in DATA {
            sum += (PauliWord::from(word), coeff);
        }
        sum
    }};
}

macro_rules! support_old {
    ($sum:expr) => {{
        let mut out: Vec<_> = $sum.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }};
}

macro_rules! support_new {
    ($sum:expr) => {{
        let mut out: Vec<_> = $sum.iter().map(|(k, c)| (k.to_string(), c)).collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }};
}

macro_rules! bench_truncate {
    ($c:expr, $name:expr, $old:expr, $new:expr) => {{
        let (old, new) = ($old, $new);
        let (mut op, mut np) = (old.clone(), new.clone());
        op.truncate();
        np.truncate();
        assert_eq!(support_old!(op), support_new!(np), "{}", $name);
        let mut group = $c.benchmark_group(concat!("pauli_sum_surface/truncate/", $name));
        group.bench_function("old", |b| {
            b.iter_batched_ref(|| old.clone(), |s| s.truncate(), BatchSize::LargeInput)
        });
        group.bench_function("new", |b| {
            b.iter_batched_ref(|| new.clone(), |s| s.truncate(), BatchSize::LargeInput)
        });
        group.finish();
    }};
}

pub fn bench(c: &mut Criterion) {
    coefficient(c);
    weight(c);
    combined(c);
    preserve(c);
}

fn coefficient(c: &mut Criterion) {
    bench_truncate!(
        c,
        "coefficient_disabled",
        fill_old!(OldCoeffSum, OldCoeff(0.0)),
        fill_new!(NewCoeffSum, NewCoeff { threshold: 0.0 })
    );
    bench_truncate!(
        c,
        "coefficient_active",
        fill_old!(OldCoeffSum, OldCoeff(1e-6)),
        fill_new!(NewCoeffSum, NewCoeff { threshold: 1e-6 })
    );
}

fn weight(c: &mut Criterion) {
    bench_truncate!(
        c,
        "max_weight_disabled",
        fill_old!(OldWeightSum, OldWeight(usize::MAX)),
        fill_new!(NewWeightSum, NewWeight(usize::MAX))
    );
    bench_truncate!(
        c,
        "max_weight_active",
        fill_old!(OldWeightSum, OldWeight(3)),
        fill_new!(NewWeightSum, NewWeight(3))
    );
}

fn combined(c: &mut Criterion) {
    bench_truncate!(
        c,
        "combined_disabled",
        fill_old!(
            OldCombinedSum,
            OldCombined(OldCoeff(0.0), OldWeight(usize::MAX))
        ),
        fill_new!(
            NewCombinedSum,
            NewCombined(NewCoeff { threshold: 0.0 }, NewWeight(usize::MAX))
        )
    );
    bench_truncate!(
        c,
        "combined_active",
        fill_old!(OldCombinedSum, OldCombined(OldCoeff(1e-6), OldWeight(3))),
        fill_new!(
            NewCombinedSum,
            NewCombined(NewCoeff { threshold: 1e-6 }, NewWeight(3))
        )
    );
}

fn preserve(c: &mut Criterion) {
    bench_truncate!(
        c,
        "preserve_empty",
        fill_old!(OldCoeffSum, OldCoeff(1e-6)),
        fill_new!(NewCoeffSum, NewCoeff { threshold: 1e-6 })
    );

    let mut old: OldCoeffSum = PauliSum::builder()
        .n_qubits(N)
        .strategy(OldCoeff(1e-6))
        .capacity(CAPACITY)
        .preserve_strings([OldKey::from("IIIIIIII"), OldKey::from("XYZIIIII")].into())
        .build();
    let mut new = NewCoeffSum::with_capacity(N, NewCoeff { threshold: 1e-6 }, CAPACITY)
        .preserving([NewKey::from("IIIIIIII"), NewKey::from("XYZIIIII")]);
    for (word, coeff) in DATA {
        old += (word, coeff);
        new += (NewKey::from(word), coeff);
    }
    bench_truncate!(c, "preserve_nonempty", old, new);
}
