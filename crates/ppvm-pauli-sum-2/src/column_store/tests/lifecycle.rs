// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use ppvm_traits_2::Indexable;

fn numbered_word(value: usize) -> PauliWord {
    let mut word = PauliWord::new(16);
    for bit in 0..9 {
        word.set_x_bit(bit, value >> bit & 1 != 0);
    }
    word
}

fn numbered_store(len: usize) -> Store {
    let mut store = Store::with_capacity(len);
    for i in 0..len {
        store.add_term(numbered_word(i), i as f64);
    }
    store
}

#[test]
fn eighteen_of_408_dead_rows_do_not_compact() {
    let mut store = numbered_store(408);
    Retain::retain(&mut store, |_, coeff| *coeff >= 18.0);

    assert_eq!(store.primary.rows(), 408);
    assert_eq!(store.primary.len(), 390);
    assert_eq!(Support::len(&store), 390);
    assert!(!store.primary.is_live(17));
    assert!(store.primary.is_live(18));
    assert_eq!(Support::get(&store, &numbered_word(17)), None);
    assert_eq!(
        store
            .primary
            .find_any(&numbered_word(17), numbered_word(17).key_hash()),
        Some(17)
    );
}

#[test]
fn compaction_runs_at_ceil_one_eighth_and_is_stable() {
    let mut store = numbered_store(408);
    Retain::retain(&mut store, |_, coeff| *coeff >= 51.0);

    assert_eq!(store.primary.rows(), 357);
    assert_eq!(store.primary.len(), 357);
    assert!(store.primary.is_dense());
    let keys: Vec<_> = Support::iter(&store).map(|(key, _)| key).collect();
    let expected: Vec<_> = (51..408).map(numbered_word).collect();
    assert_eq!(keys, expected);
    for i in 51..408 {
        assert_eq!(Support::get(&store, &numbered_word(i)), Some(i as f64));
    }
}

#[test]
fn readded_key_appends_instead_of_reviving_its_old_row() {
    let mut store = numbered_store(16);
    let key = numbered_word(3);
    Retain::retain(&mut store, |candidate, _| candidate != &key);
    assert_eq!(store.primary.rows(), 16);

    store.add_term(key, 7.0);

    assert_eq!(store.primary.rows(), 17);
    assert_eq!(store.primary.len(), 16);
    assert!(!store.primary.is_live(3));
    assert!(store.primary.is_live(16));
    assert_eq!(store.primary.key(16), key);
    assert_eq!(Support::get(&store, &key), Some(7.0));
    assert_eq!(Support::iter(&store).last(), Some((key, 7.0)));
}

#[test]
fn stale_index_entry_is_repointed_for_add_and_insert() {
    for insert in [false, true] {
        let mut store = numbered_store(16);
        let key = numbered_word(5);
        let hash = key.key_hash();
        Retain::retain(&mut store, |candidate, _| candidate != &key);
        assert_eq!(store.primary.find_any(&key, hash), Some(5));

        if insert {
            store.insert_term(key, 11.0);
        } else {
            store.add_term(key, 11.0);
        }

        assert_eq!(store.primary.find_any(&key, hash), Some(16));
        assert_eq!(store.primary.find(&key, hash), Some(16));
        assert_eq!(Support::get(&store, &key), Some(11.0));
    }
}

#[test]
fn exact_zero_stays_live_until_explicitly_removed() {
    let mut store = numbered_store(16);
    let zero = numbered_word(0);
    let removed = numbered_word(1);
    Retain::retain(&mut store, |candidate, _| candidate != &removed);

    assert_eq!(Support::get(&store, &zero), Some(0.0));
    assert!(store.primary.is_live(0));
    assert_eq!(Support::len(&store), 15);

    Accumulate::reduce(&mut store);
    assert_eq!(Support::get(&store, &zero), None);
    assert_eq!(Support::len(&store), 14);
}

#[test]
fn clone_preserves_tombstones_and_reset_clears_all_rows() {
    let mut store = numbered_store(16);
    Retain::retain(&mut store, |_, coeff| *coeff != 2.0);
    let mut cloned = store.clone();

    assert_eq!(cloned, store);
    assert_eq!(cloned.primary.rows(), 16);
    assert_eq!(cloned.primary.len(), 15);
    assert!(!cloned.primary.is_live(2));
    cloned.primary.debug_assert_valid();

    StoreAlloc::reset(&mut cloned);
    assert_eq!(Support::len(&cloned), 0);
    assert_eq!(cloned.primary.rows(), 0);
    assert!(cloned.primary.live.is_empty());
    cloned.primary.debug_assert_valid();
}

#[test]
fn iteration_pairing_and_scale_visit_live_rows_in_physical_order() {
    let mut store = numbered_store(32);
    Retain::retain(&mut store, |_, coeff| *coeff != 1.0 && *coeff != 4.0);
    assert_eq!(store.primary.rows(), 32);

    let before_dead = [store.primary.coeffs[1], store.primary.coeffs[4]];
    Scale::scale(&mut store, &2.0);
    assert_eq!(
        [store.primary.coeffs[1], store.primary.coeffs[4]],
        before_dead
    );

    let got: Vec<_> = Support::iter(&store).collect();
    let expected: Vec<_> = (0..32)
        .filter(|&i| i != 1 && i != 4)
        .map(|i| (numbered_word(i), 2.0 * i as f64))
        .collect();
    assert_eq!(got, expected);

    let expected_overlap: f64 = expected.iter().map(|(_, coeff)| coeff * coeff).sum();
    assert_eq!(Pair::overlap(&store, &store), expected_overlap);
}
