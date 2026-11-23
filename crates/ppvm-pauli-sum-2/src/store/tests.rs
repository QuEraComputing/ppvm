// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_lossy_pauli_word_2::LossyPauliWord;
use ppvm_pauli_word_2::PauliWord;
use ppvm_traits_2::{KeyBatch, TermSink, Word};

use super::*;

type Store = HashMapStore<PauliWord, f64>;

#[test]
fn accumulation_preserves_zero_terms_until_reduce() {
    let x = PauliWord::from("X");
    let z = PauliWord::from("Z");
    let mut store = Store::with_capacity(2);
    let mut terms = TermBatch::new();
    terms.push(x, 2.0);
    terms.push(x, -2.0);
    terms.push(z, 0.0);
    store.accumulate_batch(&terms);

    assert_eq!(store.len(), 2);
    assert_eq!(store.get(&x), Some(0.0));
    assert_eq!(store.get(&z), Some(0.0));

    terms.clear();
    terms.push(x, 3.0);
    store.accumulate_batch(&terms);
    store.accumulate(x, -1.0);
    store.reduce();
    assert_eq!(store.len(), 1);
    assert_eq!(store.get(&x), Some(2.0));
    assert_eq!(store.get(&z), None);
}

#[test]
fn scaling_and_retaining_do_not_implicitly_reduce() {
    let mut store = Store::with_capacity(3);
    for (word, coeff) in [("II", 1.0), ("XI", 2.0), ("XY", 3.0)] {
        store.accumulate(word.into(), coeff);
    }
    store.scale(&2.0);
    store.retain(|key, coeff| key.weight() <= 1 && *coeff >= 4.0);
    assert_eq!(store.len(), 1);
    assert_eq!(store.get(&"XI".into()), Some(4.0));

    store.scale(&0.0);
    assert_eq!(store.get(&"XI".into()), Some(0.0));
    store.reduce();
    assert!(store.is_empty());
}

#[test]
fn branching_merges_many_branches_per_term_after_all_visits() {
    let mut store = Store::with_capacity(4);
    for (word, coeff) in [("X", 1.0), ("Y", 2.0), ("Z", 3.0)] {
        store.accumulate(word.into(), coeff);
    }
    let (x, y) = (PauliWord::from("X"), PauliWord::from("Y"));
    let mut visits = 0;
    store.branch_in_place(|key, coeff| {
        visits += 1;
        *coeff *= 10.0;
        if *key == x {
            // Two branches: one new key, one colliding with an original key.
            vec![("I".into(), 0.5), ("Z".into(), 4.0)]
        } else if *key == y {
            // Cancels X's new branch; the exact zero stays stored.
            vec![("I".into(), -0.5)]
        } else {
            vec![]
        }
    });
    // Only originals are visited, and Z is scaled before X's branch lands:
    // merging early would give (3 + 4) * 10 = 70.
    assert_eq!(visits, 3);
    assert_eq!(store.len(), 4);
    for (word, coeff) in [("I", 0.0), ("X", 10.0), ("Y", 20.0), ("Z", 34.0)] {
        assert_eq!(store.get(&word.into()), Some(coeff), "{word}");
    }
    assert!(store.buffer.is_empty());
}

#[test]
fn complex_pairings_preserve_operand_order() {
    let x = PauliWord::from("X");
    let mut left = HashMapStore::<PauliWord, Complex<f64>>::with_capacity(2);
    left.accumulate(x, Complex::new(1.0, 2.0));
    left.accumulate(PauliWord::from("Z"), Complex::new(5.0, 0.0));
    let mut right = HashMapStore::with_capacity(1);
    right.accumulate(x, Complex::new(3.0, 4.0));

    assert_eq!(left.overlap(&right), Complex::new(-5.0, 10.0));
    assert_eq!(right.overlap(&left), Complex::new(-5.0, 10.0));
    assert_eq!(left.hermitian_overlap(&right), Complex::new(11.0, -2.0));
    assert_eq!(right.hermitian_overlap(&left), Complex::new(11.0, 2.0));

    let empty = HashMapStore::with_capacity(0);
    assert_eq!(left.overlap(&empty), Complex::new(0.0, 0.0));
    assert_eq!(empty.hermitian_overlap(&left), Complex::new(0.0, 0.0));
}

#[test]
fn support_iteration_and_batch_probes_agree() {
    let x = PauliWord::from("X");
    let z = PauliWord::from("Z");
    let mut store = Store::with_capacity(2);
    store.accumulate(x, 2.0);
    store.accumulate(z, 0.0);

    let mut visited = Vec::new();
    store.for_each_ref(|key, coeff| visited.push((*key, *coeff)));
    assert_eq!(visited, store.iter().collect::<Vec<_>>());

    let mut keys = KeyBatch::new();
    keys.push(z);
    keys.push(PauliWord::from("I"));
    keys.push(x);
    keys.fill_hashes();
    let mut values = [None; 3];
    store.probe_batch(&keys, &mut values);
    assert_eq!(values, [Some(0.0), None, Some(2.0)]);
}

#[test]
fn lossy_words_keep_distinct_support() {
    let mut store = HashMapStore::<LossyPauliWord, f64>::with_capacity(2);
    store.accumulate("LI".into(), 2.0);
    store.accumulate("II".into(), 3.0);
    store.accumulate("LI".into(), -1.0);
    assert_eq!(store.len(), 2);
    assert_eq!(store.get(&"LI".into()), Some(1.0));
    store.retain(|key, _| key.loss_weight() == 0);
    assert_eq!(store.len(), 1);
    assert_eq!(store.get(&"II".into()), Some(3.0));
}

#[test]
fn clone_and_clear_keep_storage_independent() {
    let x = PauliWord::from("X");
    let mut store = Store::with_capacity(8);
    store.accumulate(x, 2.0);
    let capacity = store.terms.capacity();

    let mut cloned = store.clone();
    assert_eq!(cloned.get(&x), Some(2.0));
    cloned.accumulate(x, 1.0);
    assert_eq!(store.get(&x), Some(2.0));

    store.clear();
    assert!(store.is_empty());
    assert_eq!(store.terms.capacity(), capacity);
    assert_eq!(cloned.get(&x), Some(3.0));
    store.accumulate(x, 4.0);
    assert_eq!(store.get(&x), Some(4.0));
}

#[test]
fn hash_collisions_do_not_merge_distinct_keys() {
    #[derive(Clone, PartialEq, Eq)]
    struct Key(u8);

    impl std::hash::Hash for Key {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            state.write_u64(self.key_hash());
        }
    }

    impl Indexable for Key {
        fn key_hash(&self) -> u64 {
            0
        }
    }

    let mut store = HashMapStore::with_capacity(2);
    store.accumulate(Key(0), 1.0);
    store.accumulate(Key(1), 2.0);
    let mut batch = TermBatch::new();
    batch.push(Key(0), 3.0);
    batch.push(Key(2), 5.0);
    store.accumulate_batch(&batch);
    assert_eq!(store.len(), 3);
    assert_eq!(store.get(&Key(0)), Some(4.0));
    assert_eq!(store.get(&Key(1)), Some(2.0));
    assert_eq!(store.get(&Key(2)), Some(5.0));
}

#[test]
fn partial_rekey_matches_word_conjugation() {
    use crate::PauliSum;
    use ppvm_phased_pauli_word_2::Phased;
    use ppvm_traits_2::{Clifford, CliffordBatch, Phase};

    // All 4-qubit words, so most gates leave some terms in place and move others.
    let words: Vec<(PauliWord, f64)> = (0..256)
        .map(|i| {
            let s: String = (0..4)
                .map(|q| ['I', 'X', 'Y', 'Z'][(i >> (2 * q)) & 3])
                .collect();
            (PauliWord::from(s.as_str()), i as f64 + 1.0)
        })
        .collect();
    type SumGate = fn(&mut PauliSum);
    type WordGate = fn(&mut Phased<PauliWord>);
    let gates: [(SumGate, WordGate); 3] = [
        (|s| s.h(1), |w| w.h(1)),
        (|s| s.cnot(2, 0), |w| w.cnot(2, 0)),
        (
            |s| s.s_many(&[0, 3]),
            |w| {
                w.s(0);
                w.s(3);
            },
        ),
    ];
    for (apply_sum, apply_word) in gates {
        let mut sum = PauliSum::from_terms(4, words.clone());
        apply_sum(&mut sum);
        assert_eq!(sum.len(), words.len());
        for (word, coeff) in words.as_slice() {
            let mut phased = Phased::new(*word);
            apply_word(&mut phased);
            let (key, phase) = phased.into_parts();
            let expected = if phase == Phase::Neg1 { -coeff } else { *coeff };
            assert_eq!(sum.get(&key), Some(expected));
        }
    }

    let original: PauliSum = PauliSum::from_terms(4, words.clone());
    let mut twice = original.clone();
    twice.h(2);
    twice.h(2);
    assert_eq!(twice, original);
}
