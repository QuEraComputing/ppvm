// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use ppvm_pauli_word_2::PauliWord;
use ppvm_traits_2::PauliBits;

type Store = ColumnStore<PauliWord, f64>;

mod lifecycle;
mod rotations;

fn store(terms: &[(&str, f64)]) -> Store {
    let mut store = Store::with_capacity(8);
    for (word, coeff) in terms {
        store.add_term(PauliWord::from(*word), *coeff);
    }
    store
}

#[test]
fn add_term_accumulates_and_keeps_zero() {
    let s = store(&[("XI", 1.0), ("IZ", 2.0), ("XI", 3.0)]);
    assert_eq!(Support::len(&s), 2);
    assert_eq!(Support::get(&s, &PauliWord::from("XI")), Some(4.0));

    let s = store(&[("XI", 1.0), ("XI", -1.0)]);
    assert_eq!(Support::len(&s), 1);
    assert_eq!(Support::get(&s, &PauliWord::from("XI")), Some(0.0));
}

#[test]
fn reduce_is_the_only_zero_dropper() {
    let mut store = store(&[("XI", 1.0), ("IZ", 0.0)]);
    assert_eq!(Support::len(&store), 2);
    Accumulate::reduce(&mut store);
    assert_eq!(Support::len(&store), 1);
    assert_eq!(Support::get(&store, &PauliWord::from("XI")), Some(1.0));
    assert_eq!(Support::get(&store, &PauliWord::from("IZ")), None);
}

// The `AtomicU64` hash cache is interior-mutable but excluded from `Eq`/`Hash`.
#[test]
#[allow(clippy::mutable_key_type)]
fn probe_survives_index_growth() {
    let mut store = Store::with_capacity(0);
    let words: Vec<PauliWord> = (0..512u32)
        .map(|i| {
            let mut word = PauliWord::new(16);
            for bit in 0..16 {
                if i >> (bit % 9) & 1 == 1 {
                    word.set_x_bit(bit, true);
                }
                if i >> (bit % 7) & 1 == 1 {
                    word.set_z_bit(bit, true);
                }
            }
            word
        })
        .collect();
    for (i, word) in words.iter().enumerate() {
        store.add_term(word.clone(), i as f64);
    }
    for word in &words {
        assert!(
            Support::get(&store, word).is_some(),
            "lost key after growth"
        );
    }
    let distinct: std::collections::HashSet<_> = words.iter().collect();
    assert_eq!(Support::len(&store), distinct.len());
}

#[test]
fn retain_compacts_and_keeps_probing() {
    let mut store = store(&[("XI", 1.0), ("IZ", 0.5), ("YY", 2.0), ("ZZ", 0.25)]);
    Retain::retain(&mut store, |_, coeff| *coeff >= 1.0);
    assert_eq!(Support::len(&store), 2);
    assert_eq!(Support::get(&store, &PauliWord::from("XI")), Some(1.0));
    assert_eq!(Support::get(&store, &PauliWord::from("YY")), Some(2.0));
    assert_eq!(Support::get(&store, &PauliWord::from("IZ")), None);
    store.add_term(PauliWord::from("IZ"), 7.0);
    assert_eq!(Support::get(&store, &PauliWord::from("IZ")), Some(7.0));
    assert_eq!(Support::len(&store), 3);
}

#[test]
fn rekey_in_place_preserves_coefficients() {
    let mut store = store(&[("XI", 1.0), ("IZ", 2.0)]);
    store.rekey_bijective(|key, coeff| {
        let mut new_key = PauliWord::new(2);
        new_key.set_x_bit(1, key.x_bit(0));
        new_key.set_z_bit(1, key.z_bit(0));
        new_key.set_x_bit(0, key.x_bit(1));
        new_key.set_z_bit(0, key.z_bit(1));
        (new_key, coeff)
    });
    assert_eq!(Support::len(&store), 2);
    assert_eq!(Support::get(&store, &PauliWord::from("IX")), Some(1.0));
    assert_eq!(Support::get(&store, &PauliWord::from("ZI")), Some(2.0));
}

#[test]
fn rotate_is_two_pass() {
    let mut store = store(&[("Z", 1.0), ("Y", 1.0)]);
    store.rotate_in_place(|key, coeff| {
        let branch_is_y = key.z_bit(0) && !key.x_bit(0);
        let out = if branch_is_y {
            PauliWord::from("Y")
        } else {
            PauliWord::from("Z")
        };
        *coeff *= 0.5;
        Some((out, 0.25))
    });
    assert_eq!(Support::get(&store, &PauliWord::from("Z")), Some(0.75));
    assert_eq!(Support::get(&store, &PauliWord::from("Y")), Some(0.75));
}

#[test]
fn clone_drops_workspace_and_equality_ignores_order() {
    let a = store(&[("XI", 1.0), ("IZ", 2.0)]);
    let b = store(&[("IZ", 2.0), ("XI", 1.0)]);
    assert_eq!(a, b, "insertion order is not part of the value");
    assert_eq!(a, a.clone());

    let c = store(&[("XI", 1.0), ("IZ", 2.0), ("YY", 0.0)]);
    assert_ne!(a, c);
}

#[test]
fn clone_preserves_reserved_workspace_capacity() {
    let mut original = Store::with_capacity(64);
    original.add_term(PauliWord::from("XI"), 1.0);
    original.scratch.reserve(17);
    original.batch = TermBatch::with_capacity(23);

    let cloned = original.clone();
    assert!(cloned.primary.keys.capacity() >= 64);
    assert!(cloned.primary.coeffs.capacity() >= 64);
    assert!(cloned.aux.keys.capacity() >= 64);
    assert!(cloned.aux.coeffs.capacity() >= 64);
    assert!(cloned.scratch.capacity() >= original.scratch.capacity());
    assert!(cloned.batch.capacity() >= original.batch.capacity());
    assert_eq!(cloned, original);
}
