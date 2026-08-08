// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::{BatchSize, Criterion};
use ppvm_pauli_sum::config::fxhash::Byte;
use ppvm_pauli_sum::sum::PauliSum;
use ppvm_pauli_sum_2::{HashMapStore, NoPolicy, PauliWord, Sum};
use ppvm_pauli_word::word::PauliWord as OldPauliWord;
use ppvm_traits::traits::NoStrategy;

pub const N: usize = 8;
pub const CAPACITY: usize = 1024;
const SUPPORT: usize = 192;

pub type OldKey = OldPauliWord<[u8; 8]>;
pub type NewKey = PauliWord<[u8; 8]>;
pub type OldConfig = Byte<8, f64, NoStrategy, OldKey>;
pub type OldSum = PauliSum<OldConfig>;
pub type NewSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;

pub fn terms(offset: usize, count: usize) -> Vec<(String, f64)> {
    const LETTERS: [char; 4] = ['I', 'X', 'Y', 'Z'];
    (offset..offset + count)
        .map(|i| {
            let mut value = i.wrapping_mul(4051).wrapping_add(17) & 0xffff;
            let word = (0..N)
                .map(|_| {
                    let p = LETTERS[value & 3];
                    value >>= 2;
                    p
                })
                .collect();
            let coeff = 0.25 + ((i * 29 % 97) as f64) / 31.0;
            (word, if i & 1 == 0 { coeff } else { -coeff })
        })
        .collect()
}

pub fn keyed_old(data: &[(String, f64)]) -> Vec<(OldKey, f64)> {
    data.iter()
        .map(|(w, c)| (OldKey::from(w.as_str()), *c))
        .collect()
}

pub fn keyed_new(data: &[(String, f64)]) -> Vec<(NewKey, f64)> {
    data.iter()
        .map(|(w, c)| (NewKey::from(w.as_str()), *c))
        .collect()
}

pub fn build_old(data: &[(String, f64)]) -> OldSum {
    let mut sum = OldSum::builder()
        .n_qubits(N)
        .strategy(NoStrategy)
        .capacity(CAPACITY)
        .build();
    for (key, coeff) in keyed_old(data) {
        sum += (key, coeff);
    }
    sum
}

pub fn build_new(data: &[(String, f64)]) -> NewSum {
    let mut sum = NewSum::with_capacity(N, NoPolicy, CAPACITY);
    for (key, coeff) in keyed_new(data) {
        sum += (key, coeff);
    }
    sum
}

pub fn prepared() -> (OldSum, NewSum) {
    let data = terms(0, SUPPORT);
    (build_old(&data), build_new(&data))
}

#[track_caller]
pub fn assert_pair(old: &OldSum, new: &NewSum) {
    let mut old_terms: Vec<_> = old.iter().map(|(k, c)| (k.to_string(), *c)).collect();
    let mut new_terms: Vec<_> = new.iter().map(|(k, c)| (k.to_string(), c)).collect();
    old_terms.sort_by(|a, b| a.0.cmp(&b.0));
    new_terms.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(old_terms.len(), new_terms.len(), "support size differs");
    for ((ok, oc), (nk, nc)) in old_terms.iter().zip(&new_terms) {
        assert_eq!(ok, nk, "support key differs");
        let tol = 1e-9_f64.max(oc.abs() * 1e-10);
        assert!((oc - nc).abs() <= tol, "{ok}: old={oc}, new={nc}");
    }
}

pub fn bench_mut<FO, FN>(c: &mut Criterion, name: &str, old_op: FO, new_op: FN)
where
    FO: Fn(&mut OldSum) + Copy,
    FN: Fn(&mut NewSum) + Copy,
{
    let (old_seed, new_seed) = prepared();
    assert_pair(&old_seed, &new_seed);
    let (mut old_probe, mut new_probe) = (old_seed.clone(), new_seed.clone());
    old_op(&mut old_probe);
    new_op(&mut new_probe);
    assert_pair(&old_probe, &new_probe);

    let mut group = c.benchmark_group(format!("pauli_sum_surface/{name}"));
    group.bench_function("old", |b| {
        b.iter_batched_ref(|| old_seed.clone(), old_op, BatchSize::LargeInput)
    });
    group.bench_function("new", |b| {
        b.iter_batched_ref(|| new_seed.clone(), new_op, BatchSize::LargeInput)
    });
    group.finish();
}
