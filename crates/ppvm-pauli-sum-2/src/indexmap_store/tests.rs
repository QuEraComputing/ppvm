// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_traits_2::{Accumulate, Multiply, Pair, Support, TermBatch, TermSink};

use super::IndexMapStore;
use crate::store::{
    AddTerm, ApplyProducer, BranchInPlace, InsertTerm, MultiplyInPlace, RekeyBijective,
    RotateInPlace, StoreAlloc,
};
use crate::{IndexPauliSum, NoPolicy, PauliWord, RekeyProducer};

type Key = PauliWord<[u8; 1]>;
type Store = IndexMapStore<Key, f64>;

fn key(word: &str) -> Key {
    Key::from(word)
}

fn order(store: &Store) -> Vec<String> {
    Support::iter(store).map(|(k, _)| k.to_string()).collect()
}

#[test]
fn add_replace_reduce_preserve_order_and_exact_zero() {
    let mut store = Store::with_capacity(8);
    AddTerm::add_term(&mut store, key("X"), 1.0);
    AddTerm::add_term(&mut store, key("Z"), 2.0);
    AddTerm::add_term(&mut store, key("X"), -1.0);
    InsertTerm::insert_term(&mut store, key("X"), 7.0);
    AddTerm::add_term(&mut store, key("Y"), 0.0);

    assert_eq!(order(&store), ["X", "Z", "Y"]);
    assert_eq!(Support::get(&store, &key("X")), Some(7.0));
    assert_eq!(Support::get(&store, &key("Y")), Some(0.0));

    Accumulate::reduce(&mut store);
    assert_eq!(order(&store), ["X", "Z"]);
}

#[test]
fn extend_replaces_without_moving_existing_key() {
    let mut sum: IndexPauliSum<1> = IndexPauliSum::new(3);
    sum.extend([(key("XII"), 1.0), (key("ZII"), 2.0)]);
    sum.extend([(key("XII"), 9.0), (key("YII"), 3.0)]);
    let terms: Vec<_> = sum.iter().map(|(k, c)| (k.to_string(), c)).collect();
    assert_eq!(
        terms,
        [
            ("XII".into(), 9.0),
            ("ZII".into(), 2.0),
            ("YII".into(), 3.0)
        ]
    );
}

#[test]
fn rekey_and_rotation_follow_legacy_index_order() {
    let mut store = Store::with_capacity(4);
    for (word, coeff) in [("X", 1.0), ("Z", 2.0), ("I", 3.0)] {
        AddTerm::add_term(&mut store, key(word), coeff);
    }
    RekeyBijective::rekey_bijective(&mut store, |k, c| {
        let next = match k.to_string().as_str() {
            "X" => key("Y"),
            "Z" => key("X"),
            _ => key("Z"),
        };
        (next, c)
    });
    assert_eq!(order(&store), ["Y", "X", "Z"]);

    RotateInPlace::rotate_in_place(&mut store, |k, c| {
        *c *= 10.0;
        match k.to_string().as_str() {
            "Y" => Some((key("X"), 1.0)), // collision keeps X's position
            "X" => Some((key("I"), 0.0)), // new zero appends
            _ => None,
        }
    });
    assert_eq!(order(&store), ["Y", "X", "Z", "I"]);
    assert_eq!(Support::get(&store, &key("X")), Some(21.0));
    assert_eq!(Support::get(&store, &key("I")), Some(0.0));
}

#[test]
fn high_fanout_branch_places_branch_map_before_unique_survivors() {
    let mut store = Store::with_capacity(2);
    AddTerm::add_term(&mut store, key("X"), 1.0);
    AddTerm::add_term(&mut store, key("Z"), 2.0);
    BranchInPlace::branch_in_place(&mut store, |k, _, sink| {
        if k == &key("X") {
            sink.extend([(key("I"), 3.0), (key("Y"), 4.0), (key("X"), 5.0)]);
        }
    });
    assert_eq!(order(&store), ["I", "Y", "X", "Z"]);
    assert_eq!(Support::get(&store, &key("X")), Some(6.0));
}

#[test]
fn duplicate_fanout_uses_deduplicated_branch_cardinality() {
    let mut store = Store::with_capacity(4);
    AddTerm::add_term(&mut store, key("X"), 1.0);
    AddTerm::add_term(&mut store, key("Z"), 2.0);
    BranchInPlace::branch_in_place(&mut store, |k, _, sink| {
        if k == &key("X") {
            // Raw fan-out 4 > primary len 2, but only two unique branch keys.
            sink.extend([
                (key("I"), 1.0),
                (key("I"), 2.0),
                (key("Y"), 3.0),
                (key("Y"), 4.0),
            ]);
        }
    });
    assert_eq!(order(&store), ["X", "Z", "I", "Y"]);
    assert_eq!(Support::get(&store, &key("I")), Some(3.0));
    assert_eq!(Support::get(&store, &key("Y")), Some(7.0));
}

#[test]
fn all_workspace_capacities_survive_reuse() {
    let mut store = Store::with_capacity(16);
    AddTerm::add_term(&mut store, key("X"), 1.0);
    RotateInPlace::rotate_in_place(&mut store, |_, _| Some((key("Z"), 0.0)));
    ApplyProducer::apply_producer(&mut store, RekeyProducer::new(|k: &Key, c: &f64| (*k, *c)));

    let capacities = (
        store.primary.capacity(),
        store.aux.capacity(),
        store.scratch.capacity(),
        store.batch.capacity(),
    );
    let cloned = store.clone();
    assert_eq!(
        capacities,
        (
            cloned.primary.capacity(),
            cloned.aux.capacity(),
            cloned.scratch.capacity(),
            cloned.batch.capacity(),
        )
    );
    StoreAlloc::reset(&mut store);
    assert!(store.primary.is_empty());
    assert_eq!(
        capacities,
        (
            store.primary.capacity(),
            store.aux.capacity(),
            store.scratch.capacity(),
            store.batch.capacity(),
        )
    );

    let mut batch = TermBatch::new();
    batch.push(key("Y"), 0.0);
    Accumulate::accumulate_batch(&mut store, &batch);
    assert_eq!(order(&store), ["Y"]);
}

#[test]
fn clone_and_equality_ignore_transient_workspace_contents() {
    let mut store = Store::with_capacity(8);
    AddTerm::add_term(&mut store, key("X"), 1.0);
    let cloned = store.clone();
    assert_eq!(store, cloned);
    assert_eq!(order(&store), order(&cloned));
    assert!(cloned.aux.is_empty());
    assert!(cloned.scratch.is_empty());
    assert!(cloned.batch.is_empty());
    let _: IndexPauliSum<1, f64, NoPolicy> = IndexPauliSum::new(1);
}

#[test]
fn pairing_and_products_cover_the_complete_algebra_surface() {
    let mut left = Store::with_capacity(4);
    let mut right = Store::with_capacity(4);
    for (word, coeff) in [("X", 2.0), ("Z", 3.0)] {
        AddTerm::add_term(&mut left, key(word), coeff);
    }
    for (word, coeff) in [("Z", 5.0), ("X", 7.0)] {
        AddTerm::add_term(&mut right, key(word), coeff);
    }
    assert_eq!(Pair::overlap(&left, &right), 29.0);
    assert_eq!(Pair::hermitian_overlap(&left, &right), 29.0);

    type ComplexStore = IndexMapStore<Key, Complex<f64>>;
    let mut a = ComplexStore::with_capacity(4);
    let mut b = ComplexStore::with_capacity(4);
    AddTerm::add_term(&mut a, key("X"), Complex::new(1.0, 0.0));
    AddTerm::add_term(&mut a, key("Z"), Complex::new(2.0, 0.0));
    AddTerm::add_term(&mut b, key("I"), Complex::new(3.0, 0.0));
    AddTerm::add_term(&mut b, key("X"), Complex::new(4.0, 0.0));
    let mut product = ComplexStore::with_capacity(8);
    Multiply::multiply_into(&a, &b, &mut product);
    assert_eq!(
        Support::iter(&product)
            .map(|(k, _)| k.to_string())
            .collect::<Vec<_>>(),
        ["X", "I", "Z", "Y"]
    );
    MultiplyInPlace::multiply_in_place(&mut a, &b);
    assert_eq!(a, product);
    assert_eq!(
        Support::iter(&a).collect::<Vec<_>>(),
        Support::iter(&product).collect::<Vec<_>>()
    );
}
