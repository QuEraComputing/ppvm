// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential coverage for the lossy `PauliSum` gate family. The old crate is
//! the behavioral oracle; coefficients and zero-valued support entries are
//! compared bit-exactly.

use ppvm_lossy_pauli_word_2::LossyPauliWord as NewWord;
use ppvm_pauli_sum::config;
use ppvm_pauli_sum::strategy::MaxLossWeight as OldMaxLossWeight;
use ppvm_pauli_sum::sum::PauliSum as OldSumT;
use ppvm_pauli_sum_2::{LossyPauliSum as NewSumT, MaxLossWeight, NoPolicy, Sum};
use ppvm_pauli_word::loss::LossyPauliWord as OldWord;
use ppvm_traits::traits::{
    Clifford as OldClifford, CorrelatedLossChannel as OldCorrelatedLossChannel,
    LossChannel as OldLossChannel, NoStrategy, ResetLossChannel as OldResetLossChannel,
    RotationOne as OldRotationOne, Trace as OldTrace,
};
use ppvm_traits_2::{
    Clifford as NewClifford, CorrelatedLossChannel as NewCorrelatedLossChannel,
    LossChannel as NewLossChannel, ResetLossChannel as NewResetLossChannel,
    RotationOne as NewRotationOne, Trace as NewTrace,
};

type OldConfig = config::fxhash::Byte<1, f64, NoStrategy, OldWord<[u8; 1], fxhash::FxBuildHasher>>;
type OldSum = OldSumT<OldConfig>;
type NewSum = NewSumT<f64, NoPolicy>;

fn old_sum(n: usize, terms: &[(&str, f64)]) -> OldSum {
    let mut sum = OldSum::builder().n_qubits(n).build();
    for &(word, coeff) in terms {
        sum += (word, coeff);
    }
    sum
}

fn new_sum(n: usize, terms: &[(&str, f64)]) -> NewSum {
    Sum::from_terms(
        n,
        terms
            .iter()
            .map(|&(word, coeff)| (NewWord::from(word), coeff)),
    )
}

fn old_support(sum: &OldSum) -> Vec<(String, u64)> {
    let mut out: Vec<_> = sum
        .data()
        .iter()
        .map(|(word, coeff)| (word.to_string(), coeff.to_bits()))
        .collect();
    out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out
}

fn new_support(sum: &NewSum) -> Vec<(String, u64)> {
    let mut out: Vec<_> = sum
        .iter()
        .map(|(word, coeff)| (word.to_string(), coeff.to_bits()))
        .collect();
    out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    out
}

fn assert_same(old: &OldSum, new: &NewSum) {
    assert_eq!(old_support(old), new_support(new));
}

#[test]
fn reset_and_single_loss_match_exactly_including_zero_entries() {
    for word in ["I", "X", "Y", "Z", "L"] {
        let mut old = old_sum(1, &[(word, 1.0)]);
        let mut new = new_sum(1, &[(word, 1.0)]);
        old.reset_loss_channel(0);
        new.reset_loss_channel(0);
        assert_same(&old, &new);

        old.loss_channel(0, 0.2);
        new.loss_channel(0, 0.2, &mut ppvm_conformance_2::analytic_rng());
        assert_same(&old, &new);
    }

    let mut old = old_sum(1, &[("L", 1.0)]);
    let mut new = new_sum(1, &[("L", 1.0)]);
    old.reset_loss_channel(0);
    new.reset_loss_channel(0);
    assert_same(&old, &new);
    assert_eq!(new_support(&new), vec![("L".into(), 0.0f64.to_bits())]);
}

#[test]
fn correlated_loss_matches_all_four_loss_arms() {
    let terms = [
        ("II", 1.0),
        ("LI", 2.0),
        ("IL", 3.0),
        ("LL", 4.0),
        ("XZ", -0.5),
    ];
    let probabilities = [0.07, 0.11, 0.19];
    let mut old = old_sum(2, &terms);
    let mut new = new_sum(2, &terms);
    old.correlated_loss_channel(0, 1, probabilities);
    new.correlated_loss_channel(0, 1, probabilities, &mut ppvm_conformance_2::analytic_rng());
    assert_same(&old, &new);
}

#[test]
fn loss_interleaved_ghz_workload_matches() {
    let mut old = old_sum(2, &[("ZZ", 1.0)]);
    let mut new = new_sum(2, &[("ZZ", 1.0)]);

    for q in 0..2 {
        old.reset_loss_channel(q);
        new.reset_loss_channel(q);
        old.loss_channel(q, 0.1);
        new.loss_channel(q, 0.1, &mut ppvm_conformance_2::analytic_rng());
    }
    old.cnot(0, 1);
    new.cnot(0, 1);
    old.loss_channel(0, 0.2);
    new.loss_channel(0, 0.2, &mut ppvm_conformance_2::analytic_rng());
    old.h(0);
    new.h(0);
    old.rx(1, 0.37);
    new.rx(1, 0.37);

    assert_same(&old, &new);
    let old_pattern: ppvm_pauli_word::pattern::PauliPattern = "Z?*".into();
    let new_pattern = ppvm_pauli_sum_2::PauliPattern::zero_state();
    let old_value: f64 = old.trace(&old_pattern);
    let new_value: f64 = new.trace(&new_pattern);
    assert_eq!(old_value.to_bits(), new_value.to_bits());
}

#[test]
fn max_loss_weight_keeps_the_same_support() {
    type OldLossConfig =
        config::fxhash::Byte<1, f64, OldMaxLossWeight, OldWord<[u8; 1], fxhash::FxBuildHasher>>;

    let mut old: OldSumT<OldLossConfig> = OldSumT::builder()
        .n_qubits(3)
        .strategy(OldMaxLossWeight(1))
        .build();
    let mut new: NewSumT<f64, MaxLossWeight> = Sum::with_policy(3, MaxLossWeight(1));
    for (word, coeff) in [("III", 1.0), ("LII", 2.0), ("LLI", 3.0), ("LLL", 4.0)] {
        old += (word, coeff);
        new += (NewWord::from(word), coeff);
    }
    old.truncate();
    new.truncate();

    let mut old_keys: Vec<_> = old.data().keys().map(ToString::to_string).collect();
    let mut new_keys: Vec<_> = new.iter().map(|(word, _)| word.to_string()).collect();
    old_keys.sort_unstable();
    new_keys.sort_unstable();
    assert_eq!(old_keys, new_keys);
}
