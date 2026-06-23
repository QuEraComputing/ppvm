// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Integration test: verify that the full MSD circuit using batch methods
//! produces identical measurement outcomes to the naive individual-gate version.

use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<Byte8F64<2>, u128>;

fn encode_naive(tab: &mut Tab, q: &[usize]) {
    for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [7, 16] {
        tab.sqrt_y_dag(q[i]);
    }
    for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [4, 10, 14, 16] {
        tab.sqrt_y_dag(q[i]);
    }
    for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
        tab.cz(q[i], q[j]);
    }
    for i in [3, 6, 9, 10, 12, 13] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
        tab.cz(q[i], q[j]);
    }
    for i in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 14] {
        tab.sqrt_y(q[i]);
    }
    for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
        tab.cz(q[i], q[j]);
    }
    for i in [0, 2, 5, 6, 8, 10, 12] {
        tab.sqrt_y_dag(q[i]);
    }
}

fn encode_batch(tab: &mut Tab, q: &[usize]) {
    tab.sqrt_y_many(&[
        q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[8], q[9], q[10], q[11], q[12], q[13], q[14],
        q[15], q[16],
    ]);
    for [i, j] in [[1, 3], [7, 10], [12, 14], [13, 16]] {
        tab.cz(q[i], q[j]);
    }
    tab.sqrt_y_dag_many(&[q[7], q[16]]);
    for [i, j] in [[4, 7], [8, 10], [11, 14], [15, 16]] {
        tab.cz(q[i], q[j]);
    }
    tab.sqrt_y_dag_many(&[q[4], q[10], q[14], q[16]]);
    for [i, j] in [[2, 4], [6, 8], [7, 9], [10, 13], [14, 16]] {
        tab.cz(q[i], q[j]);
    }
    tab.sqrt_y_many(&[q[3], q[6], q[9], q[10], q[12], q[13]]);
    for [i, j] in [[0, 2], [3, 6], [5, 8], [10, 12], [11, 13]] {
        tab.cz(q[i], q[j]);
    }
    tab.sqrt_y_many(&[
        q[1], q[2], q[3], q[4], q[6], q[7], q[8], q[9], q[11], q[12], q[14],
    ]);
    for [i, j] in [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [12, 15]] {
        tab.cz(q[i], q[j]);
    }
    tab.sqrt_y_dag_many(&[q[0], q[2], q[5], q[6], q[8], q[10], q[12]]);
}

fn run_msd_naive(seed: u64) -> String {
    let n = 85;
    let mut tab: Tab = GeneralizedTableau::new_with_seed(n, 1e-10, seed);
    let qa: Vec<usize> = (0..n).collect();
    let ql: Vec<&[usize]> = qa.chunks_exact(17).collect();

    for q in &ql {
        tab.h(q[7]);
        tab.t(q[7]);
        encode_naive(&mut tab, q);
    }
    for i in [0, 1, 4] {
        for q in ql[i] {
            tab.sqrt_x(*q);
        }
    }
    for (c, t) in ql[0].iter().zip(ql[1]) {
        tab.cz(*c, *t);
    }
    for (c, t) in ql[2].iter().zip(ql[3]) {
        tab.cz(*c, *t);
    }
    for q in ql[0] {
        tab.sqrt_y(*q);
    }
    for q in ql[3] {
        tab.sqrt_y(*q);
    }
    for (c, t) in ql[0].iter().zip(ql[2]) {
        tab.cz(*c, *t);
    }
    for (c, t) in ql[3].iter().zip(ql[4]) {
        tab.cz(*c, *t);
    }
    for q in ql[0] {
        tab.sqrt_x_dag(*q);
    }
    for (c, t) in ql[0].iter().zip(ql[4]) {
        tab.cz(*c, *t);
    }
    for (c, t) in ql[1].iter().zip(ql[3]) {
        tab.cz(*c, *t);
    }
    for block in ql.iter().take(5) {
        for q in *block {
            tab.sqrt_x_dag(*q);
        }
    }
    (0..n)
        .map(|i| if tab.measure(i).unwrap() { '1' } else { '0' })
        .collect()
}

fn run_msd_batch(seed: u64) -> String {
    let n = 85;
    let mut tab: Tab = GeneralizedTableau::new_with_seed(n, 1e-10, seed);
    let qa: Vec<usize> = (0..n).collect();
    let ql: Vec<&[usize]> = qa.chunks_exact(17).collect();

    for q in &ql {
        tab.h(q[7]);
        tab.t(q[7]);
        encode_batch(&mut tab, q);
    }
    tab.sqrt_x_many(ql[0]);
    tab.sqrt_x_many(ql[1]);
    tab.sqrt_x_many(ql[4]);
    tab.cz_block_pairs(0, 17, 17);
    tab.cz_block_pairs(34, 17, 13);
    tab.cz_block_pairs_cross_word(0, 47, 1, 0, 4);
    tab.sqrt_y_many(ql[0]);
    tab.sqrt_y_many(ql[3]);
    tab.cz_block_pairs(0, 34, 17);
    tab.cz_block_pairs_cross_word(0, 51, 1, 4, 13);
    tab.cz_block_pairs(64, 17, 4);
    tab.sqrt_x_dag_many(ql[0]);
    tab.cz_block_pairs_cross_word(0, 0, 1, 4, 17);
    tab.cz_block_pairs(17, 34, 13);
    tab.cz_block_pairs_cross_word(0, 30, 1, 0, 4);
    for block in ql.iter().take(5) {
        tab.sqrt_x_dag_many(block);
    }
    (0..n)
        .map(|i| if tab.measure(i).unwrap() { '1' } else { '0' })
        .collect()
}

#[test]
fn test_msd_batch_matches_naive() {
    for seed in 0..50 {
        let naive = run_msd_naive(seed);
        let batch = run_msd_batch(seed);
        assert_eq!(naive, batch, "Mismatch at seed={seed}");
    }
}
