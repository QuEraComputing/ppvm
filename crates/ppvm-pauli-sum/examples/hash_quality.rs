// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Hypothesis check for the trotter_storage_cliff: is `PauliWord<[u8;8]>`'s
//! cached hash poorly distributed at high fill, causing the
//! `hashbrown::insert_unique` blowup observed in flamegraphs?
//!
//! After running the full Trotter circuit at N=56 we have ~83k distinct
//! PauliWords in the map. We pull their cached hashes and check:
//!   * raw u64 collisions (should be ~0 for any good 64-bit hash)
//!   * bucket collisions for hashbrown's typical capacity (the value that
//!     actually controls insert_unique probe length)
//!   * compare both for `[u8;8]` and `[u8;16]` storage

use std::collections::HashMap;
use std::hash::{BuildHasher, Hash};

use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;
use ppvm_pauli_word::word::PauliWord;

const H: f64 = 1.0;
const DT: f64 = 0.1 / H;
const TIME: f64 = 2.0 / H;
const J: f64 = 1.0 * H;
const MIN_ABS_COEFF: f64 = 1e-6;
const NOISE: [f64; 3] = [1e-4 / 4.0; 3];

fn build_and_propagate<C, T>(n: usize) -> PauliSum<C>
where
    C: Config<Coeff = f64, Strategy = CoefficientThreshold, PauliWordType = T>,
    T: PauliWordTrait,
    for<'a> &'a str: Into<T>,
{
    let strat = CoefficientThreshold(MIN_ABS_COEFF);
    let mut state: PauliSum<C> = PauliSum::builder()
        .n_qubits(n)
        .strategy(strat)
        .capacity(n.pow(2))
        .build();
    for i in 0..n {
        let term: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
        state += (term.as_str(), 1.0);
    }
    let steps = (TIME / DT) as usize;
    let theta_zz = DT * J;
    let theta_x = DT * H;
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, NOISE);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, NOISE);
            state.truncate();
            state.pauli_error(i, NOISE);
            state.truncate();
            state.rzz(i, i + 1, theta_zz);
            state.truncate();
        }
    }
    state
}

fn analyze<S, H>(label: &str, keys: impl Iterator<Item = S>)
where
    S: Hash,
    H: BuildHasher + Default,
{
    let hasher = H::default();
    let hashes: Vec<u64> = keys.map(|k| hasher.hash_one(&k)).collect();

    let n = hashes.len();
    let mut counts: HashMap<u64, u32> = HashMap::new();
    for h in &hashes {
        *counts.entry(*h).or_insert(0) += 1;
    }
    let unique = counts.len();
    let raw_collisions: usize = counts.values().map(|c| (*c as usize) - 1).sum();

    // Simulate hashbrown bucket distribution at a typical load factor.
    // hashbrown uses power-of-two table sizes; for ~83k entries at load 0.875
    // we'd be in a 131072-bucket table.
    let bucket_bits = 17;
    let mask = (1u64 << bucket_bits) - 1;
    let mut bucket_counts: HashMap<u64, u32> = HashMap::new();
    for h in &hashes {
        *bucket_counts.entry(h & mask).or_insert(0) += 1;
    }
    let max_bucket = bucket_counts.values().max().copied().unwrap_or(0);
    let mean_load = n as f64 / (1u64 << bucket_bits) as f64;
    let variance: f64 = {
        let m = mean_load;
        let v: f64 = bucket_counts
            .values()
            .map(|c| (*c as f64 - m).powi(2))
            .sum();
        let unfilled = (1u64 << bucket_bits) as usize - bucket_counts.len();
        let zero_term = unfilled as f64 * m * m;
        (v + zero_term) / (1u64 << bucket_bits) as f64
    };

    println!(
        "{label:>26}: n={n} unique_u64={unique} raw_colls={raw_collisions}  \
         buckets@2^{bucket_bits}: load={mean_load:.3} max_bucket={max_bucket} var={variance:.3}"
    );
}

fn main() {
    let n = 56;
    println!("Propagating N={n} trotter circuit ... (this is the slow tier first)\n");

    type Cfg8 = config::indexmap::ByteFxHashF64<8, CoefficientThreshold>;
    type Cfg16 = config::indexmap::ByteFxHashF64<16, CoefficientThreshold>;

    let s8 = build_and_propagate::<Cfg8, _>(n);
    let s16 = build_and_propagate::<Cfg16, _>(n);
    println!("|s8|={} |s16|={}\n", s8.len(), s16.len());

    analyze::<&PauliWord<[u8; 8]>, fxhash::FxBuildHasher>(
        "ByteFxHashF64<8> keys",
        s8.data().keys(),
    );
    analyze::<&PauliWord<[u8; 16]>, fxhash::FxBuildHasher>(
        "ByteFxHashF64<16> keys",
        s16.data().keys(),
    );
}
