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
fn rekey_reuses_both_maps_and_clone_omits_scratch() {
    use ppvm_phased_pauli_word_2::Phased;
    use ppvm_traits_2::{Clifford, Phase};

    let mut store = Store::with_capacity(16);
    store.accumulate("X".into(), 2.0);
    store.accumulate("Y".into(), 0.0);
    let original_capacity = store.terms.capacity();
    let conjugate_h = |key, coeff: f64| {
        let mut word = Phased::new(key);
        word.h(0);
        let (key, phase) = word.into_parts();
        (key, if phase == Phase::Neg1 { -coeff } else { coeff })
    };
    store.rekey_bijective(conjugate_h);
    let rebuilt_capacity = store.terms.capacity();
    assert_eq!(store.scratch.capacity(), original_capacity);
    assert!(store.scratch.is_empty());
    assert_eq!(store.get(&"Z".into()), Some(2.0));
    assert_eq!(store.get(&"Y".into()), Some(0.0));

    let cloned = store.clone();
    assert_eq!(cloned.scratch.capacity(), 0);
    for _ in 0..3 {
        store.rekey_bijective(conjugate_h);
        assert_eq!(store.terms.capacity(), original_capacity);
        assert_eq!(store.scratch.capacity(), rebuilt_capacity);
        assert_eq!(store.get(&"X".into()), Some(2.0));
        store.rekey_bijective(conjugate_h);
        assert_eq!(store.terms.capacity(), rebuilt_capacity);
        assert_eq!(store.scratch.capacity(), original_capacity);
    }
    store.clear();
    assert_eq!(store.terms.capacity(), rebuilt_capacity);
    assert_eq!(store.scratch.capacity(), original_capacity);
    assert!(store.is_empty());
    assert_eq!(cloned.get(&"Z".into()), Some(2.0));
}

#[test]
fn sign_flips_do_not_allocate_a_second_map() {
    let mut store = Store::with_capacity(16);
    store.accumulate("X".into(), 2.0);
    store.accumulate("Z".into(), 3.0);
    let capacity = store.terms.capacity();
    store.flip_sign_by_key(|key| *key == PauliWord::from("X"));
    assert_eq!(store.get(&"X".into()), Some(-2.0));
    assert_eq!(store.get(&"Z".into()), Some(3.0));
    assert_eq!(store.terms.capacity(), capacity);
    assert_eq!(store.scratch.capacity(), 0);
}
