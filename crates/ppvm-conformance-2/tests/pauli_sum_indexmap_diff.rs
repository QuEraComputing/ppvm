// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Order-sensitive differential coverage for the Phase-7 ordered sum backend.

use std::collections::HashSet;

use ppvm_lossy_pauli_word_2::LossyPauliWord as NewLossyWord;
use ppvm_pauli_sum::config::indexmap::{ByteFxHash, ByteFxHashF64};
use ppvm_pauli_sum::strategy::CoefficientThreshold as OldThreshold;
use ppvm_pauli_sum::sum::PauliSum as OldSumT;
use ppvm_pauli_sum_2::{
    CoefficientThreshold, IndexMapStore, IndexPauliSum, NoPolicy, PauliWord, Sum,
};
use ppvm_pauli_word::loss::LossyPauliWord as OldLossyWord;
use ppvm_traits::traits::{
    Clifford as OldClifford, CorrelatedLossChannel as OldCorrelatedLoss, LossChannel as OldLoss,
    NoStrategy, RotationOne as OldRotation,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CorrelatedLossChannel as NewCorrelatedLoss, LossChannel as NewLoss,
    RotationOne as NewRotation,
};

type Old = OldSumT<ByteFxHashF64<8>>;
type New = IndexPauliSum<8>;

fn old_sum(terms: &[(&str, f64)]) -> Old {
    let n = terms.first().map_or(1, |(word, _)| word.len());
    let mut sum = Old::builder().n_qubits(n).build();
    for &(word, coeff) in terms {
        sum += (word, coeff);
    }
    sum
}

fn new_sum(terms: &[(&str, f64)]) -> New {
    let n = terms.first().map_or(1, |(word, _)| word.len());
    Sum::from_terms(
        n,
        terms
            .iter()
            .map(|&(word, coeff)| (PauliWord::from(word), coeff)),
    )
}

fn old_terms<C: ppvm_traits::config::Config>(sum: &OldSumT<C>) -> Vec<(String, f64)>
where
    for<'a> C::Map: ppvm_traits::traits::ACMapIter<'a, Item = (&'a C::PauliWordType, &'a f64)>,
    C::Coeff: Into<f64>,
    C::PauliWordType: ToString,
{
    sum.iter()
        .map(|(word, coeff)| (word.to_string(), *coeff))
        .collect()
}

fn new_terms<S, P>(sum: &Sum<S, P>) -> Vec<(String, f64)>
where
    S: ppvm_traits_2::Accumulate<Coeff = f64>,
    S::Key: ppvm_traits_2::Word + ppvm_traits_2::Indexable + ToString,
    P: ppvm_pauli_sum_2::Policy<S::Key, f64>,
{
    sum.iter()
        .map(|(word, coeff)| (word.to_string(), coeff))
        .collect()
}

fn assert_ordered_same(old: &Old, new: &New) {
    assert_eq!(old_terms(old), new_terms(new));
}

#[test]
fn terms_extend_zero_and_equal_weight_display_match() {
    let mut old = old_sum(&[("XII", 1.0), ("IZI", 2.0), ("IIZ", 3.0)]);
    let mut new = new_sum(&[("XII", 1.0), ("IZI", 2.0), ("IIZ", 3.0)]);
    old.extend([
        (ppvm_pauli_word::word::PauliWord::from("IZI"), 9.0),
        (ppvm_pauli_word::word::PauliWord::from("YII"), 0.0),
    ]);
    new.extend([(PauliWord::from("IZI"), 9.0), (PauliWord::from("YII"), 0.0)]);
    assert_ordered_same(&old, &new);
    assert_eq!(old.to_string(), new.to_string());
    assert_eq!(
        new_terms(&new),
        [
            ("XII".into(), 1.0),
            ("IZI".into(), 9.0),
            ("IIZ".into(), 3.0),
            ("YII".into(), 0.0)
        ]
    );
}

#[test]
fn clifford_rekeys_and_rotation_collisions_match_order() {
    let terms = [("XI", 1.0), ("IZ", 2.0), ("ZY", 3.0), ("YY", -4.0)];
    let mut old = old_sum(&terms);
    let mut new = new_sum(&terms);
    old.h(0);
    new.h(0);
    old.cnot(0, 1);
    new.cnot(0, 1);
    assert_ordered_same(&old, &new);

    old.rx(0, 0.37);
    new.rx(0, 0.37);
    assert_ordered_same(&old, &new);

    let mut old = old_sum(&[("Z", 1.0), ("Y", 2.0)]);
    let mut new = new_sum(&[("Z", 1.0), ("Y", 2.0)]);
    old.rx(0, std::f64::consts::FRAC_PI_2);
    new.rx(0, std::f64::consts::FRAC_PI_2);
    assert_ordered_same(&old, &new);
}

#[test]
fn truncation_and_multiple_preserve_restore_match_order() {
    type OldT = OldSumT<ByteFxHashF64<8, OldThreshold>>;
    type NewT = Sum<IndexMapStore<PauliWord<[u8; 8]>, f64>, CoefficientThreshold>;

    let keep_old = [
        ppvm_pauli_word::word::PauliWord::from("IZ"),
        ppvm_pauli_word::word::PauliWord::from("YI"),
    ];
    let mut old: OldT = OldT::builder()
        .n_qubits(2)
        .strategy(OldThreshold(0.5))
        .preserve_strings(HashSet::from(keep_old))
        .build();
    let mut new: NewT = NewT::with_policy(2, CoefficientThreshold { threshold: 0.5 })
        .preserving([PauliWord::from("IZ"), PauliWord::from("YI")]);
    for (word, coeff) in [("XI", 1.0), ("IZ", 0.1), ("YI", 0.2), ("ZI", 2.0)] {
        old += (word, coeff);
        new += (PauliWord::from(word), coeff);
    }
    old.truncate();
    new.truncate();
    assert_eq!(old_terms(&old), new_terms(&new));
    assert_eq!(
        new_terms(&new),
        [
            ("XI".into(), 1.0),
            ("ZI".into(), 2.0),
            ("IZ".into(), 0.1),
            ("YI".into(), 0.2)
        ]
    );
}

type OldLossConfig = ByteFxHash<1, f64, NoStrategy, OldLossyWord<[u8; 1], fxhash::FxBuildHasher>>;
type OldLossSum = OldSumT<OldLossConfig>;
type NewLossSum = Sum<IndexMapStore<NewLossyWord<[u8; 1]>, f64>, NoPolicy>;

#[test]
fn loss_and_correlated_loss_match_observable_order() {
    let terms = [("LL", 1.0), ("LI", 2.0), ("IL", 3.0), ("XZ", 4.0)];
    let mut old: OldLossSum = OldLossSum::builder().n_qubits(2).build();
    let mut new: NewLossSum = NewLossSum::new(2);
    for (word, coeff) in terms {
        old += (word, coeff);
        new += (NewLossyWord::from(word), coeff);
    }
    old.loss_channel(0, 0.2);
    new.loss_channel(0, 0.2, &mut ppvm_conformance_2::analytic_rng());
    old.correlated_loss_channel(0, 1, [0.07, 0.11, 0.19]);
    new.correlated_loss_channel(
        0,
        1,
        [0.07, 0.11, 0.19],
        &mut ppvm_conformance_2::analytic_rng(),
    );
    assert_eq!(old_terms(&old), new_terms(&new));
}

#[test]
fn wide_fixed_storage_keeps_order_through_gates() {
    type OldWide = OldSumT<ByteFxHashF64<32>>;
    type NewWide = IndexPauliSum<32>;
    let words = [
        format!("X{}", "I".repeat(128)),
        format!("{}Z", "I".repeat(128)),
        format!("{}Y{}", "I".repeat(64), "I".repeat(64)),
    ];
    let mut old: OldWide = OldWide::builder().n_qubits(129).capacity(8).build();
    let mut new: NewWide = NewWide::new(129);
    for (index, word) in words.iter().enumerate() {
        old += (word.as_str(), index as f64 + 1.0);
        new += (PauliWord::from(word.as_str()), index as f64 + 1.0);
    }
    old.h(0);
    new.h(0);
    old.cnot(0, 128);
    new.cnot(0, 128);
    assert_eq!(old_terms(&old), new_terms(&new));
}
