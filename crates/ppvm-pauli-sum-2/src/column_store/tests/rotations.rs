// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

fn rx_word(base: usize, partner: bool) -> PauliWord {
    let mut word = PauliWord::new(16);
    word.set_x_bit(0, partner);
    word.set_z_bit(0, true);
    for bit in 0..9 {
        word.set_x_bit(bit + 1, base >> bit & 1 != 0);
    }
    word
}

fn dead_sentinel() -> PauliWord {
    let mut word = PauliWord::new(16);
    for bit in 0..16 {
        word.set_z_bit(bit, true);
    }
    word
}

fn matched_rx_stores(pairs: usize, closed: bool) -> (Store, Store, PauliWord) {
    let mut dense = Store::with_capacity(2 * pairs);
    for base in 0..pairs {
        dense.add_term(rx_word(base, false), (2 * base + 1) as f64);
        if closed {
            dense.add_term(rx_word(base, true), (2 * base + 2) as f64);
        }
    }

    let dead = dead_sentinel();
    let mut sparse = dense.clone();
    sparse.add_term(dead.clone(), -17.0);
    Retain::retain(&mut sparse, |key, _| key != &dead);
    assert!(dense.primary.is_dense());
    assert!(!sparse.primary.is_dense());
    assert!(sparse.primary.sparse_cache_dirty);
    assert!(sparse.primary.sparse_rows.is_empty());
    assert!(sparse.primary.live_runs.is_empty());
    assert_eq!(dense, sparse);
    (dense, sparse, dead)
}

#[test]
fn dense_and_sparse_rx_match_on_large_closed_support() {
    let (mut dense, mut sparse, dead) = matched_rx_stores(300, true);
    let dead_row = sparse.primary.rows() - 1;

    RotateInPlace::rotate_x(&mut dense, 0, 0.1, 0.9);
    RotateInPlace::rotate_x(&mut sparse, 0, 0.1, 0.9);

    assert_eq!(dense, sparse);
    assert!(dense.primary.is_dense());
    assert!(!sparse.primary.sparse_cache_dirty);
    assert_eq!(sparse.primary.sparse_rows.len(), sparse.primary.len());
    assert!(!sparse.primary.is_live(dead_row));
    assert_eq!(sparse.primary.key(dead_row), dead);
    assert_eq!(sparse.primary.coeffs[dead_row], -17.0);
}

#[test]
fn dense_and_sparse_rx_match_when_support_grows() {
    let (mut dense, mut sparse, dead) = matched_rx_stores(600, false);
    let dense_before = dense.primary.rows();
    let sparse_dead_row = sparse.primary.rows() - 1;

    RotateInPlace::rotate_x(&mut dense, 0, 0.1, 0.9);
    RotateInPlace::rotate_x(&mut sparse, 0, 0.1, 0.9);

    assert_eq!(dense, sparse);
    assert_eq!(dense.primary.rows(), dense_before * 2);
    assert!(dense.primary.is_dense());
    assert!(!sparse.primary.sparse_cache_dirty);
    assert_eq!(sparse.primary.sparse_rows.len(), sparse.primary.len());
    assert!(!sparse.primary.is_live(sparse_dead_row));
    assert_eq!(sparse.primary.key(sparse_dead_row), dead);
    assert_eq!(sparse.primary.coeffs[sparse_dead_row], -17.0);
}

#[test]
fn no_op_retain_does_not_rebuild_a_dirty_sparse_cache() {
    let (_, mut sparse, _) = matched_rx_stores(32, true);
    let rows = sparse.primary.sparse_rows.clone();
    let runs = sparse.primary.live_runs.clone();

    Retain::retain(&mut sparse, |_, _| true);

    assert!(sparse.primary.sparse_cache_dirty);
    assert_eq!(sparse.primary.sparse_rows, rows);
    assert_eq!(sparse.primary.live_runs, runs);
}

#[test]
fn append_and_readd_leave_a_dirty_cache_for_lazy_rebuild() {
    let (_, mut sparse, dead) = matched_rx_stores(32, true);
    let appended = rx_word(100, false);
    let old_rows = sparse.primary.rows();

    sparse.add_term(appended.clone(), 23.0);
    sparse.add_term(dead.clone(), 29.0);

    assert!(sparse.primary.sparse_cache_dirty);
    assert!(sparse.primary.sparse_rows.is_empty());
    assert!(sparse.primary.live_runs.is_empty());
    assert_eq!(sparse.primary.rows(), old_rows + 2);
    assert_eq!(Support::get(&sparse, &appended), Some(23.0));
    assert_eq!(Support::get(&sparse, &dead), Some(29.0));

    sparse.primary.ensure_sparse_cache();
    assert!(!sparse.primary.sparse_cache_dirty);
    assert_eq!(sparse.primary.sparse_rows.len(), sparse.primary.len());
    assert_eq!(
        sparse.primary.sparse_rows[sparse.primary.len() - 2..],
        [old_rows as u32, old_rows as u32 + 1]
    );
}
