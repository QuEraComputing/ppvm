// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Cross-check Heisenberg channels against `ppvm-tableau-2` trajectories.
//! Each tableau run prepares a Z-basis (or lost) state, applies the channel
//! forward, and measures; the sum applies the channel backward to an
//! observable and evaluates it on the prepared state.

use approx::assert_abs_diff_eq;
use ppvm_pauli_sum_2::LossyPauliSum;
use ppvm_tableau_2::{ActiveDormant, DenseState, GeneralizedTableau};
use ppvm_traits_2::{
    Clifford, Depolarizing, Depolarizing2, LossChannel, Measure, PauliError, TwoQubitPauliError,
};

type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;
type ChannelPair = (fn(&mut Sim), fn(&mut LossyPauliSum));

const SHOTS: u64 = 8192;
const OBSERVABLES: [&str; 9] = ["II", "IZ", "IL", "ZI", "ZZ", "ZL", "LI", "LZ", "LL"];

/// Value of a word on a Z-basis state; `I` and `Z` vanish on lost sites.
fn basis_value(word: &str, state: &str) -> f64 {
    word.chars()
        .zip(state.chars())
        .map(|(p, s)| match (p, s) {
            ('L', 'L') => 1.0,
            (_, 'L') | ('L', _) | ('X' | 'Y', _) => 0.0,
            ('Z', '1') => -1.0,
            _ => 1.0,
        })
        .product()
}

/// Distinct two-qubit rates, so a misordered table changes the marginals.
fn distinct_rates() -> [f64; 15] {
    std::array::from_fn(|k| 0.008 * (k + 1) as f64)
}

fn assert_matches_tableau(initials: &[&str], channels: &[ChannelPair]) {
    for initial in initials {
        for &(forward, backward) in channels {
            let mut averages = [0.0; OBSERVABLES.len()];
            for seed in 0..SHOTS {
                let mut sim = Sim::new_with_seed(2, seed);
                for (q, state) in initial.chars().enumerate() {
                    match state {
                        '1' => sim.x(q),
                        'L' => sim.loss_channel(q, 1.0),
                        _ => {}
                    }
                }
                forward(&mut sim);
                let output: String = (0..2)
                    .map(|q| match sim.measure(q) {
                        None => 'L',
                        Some(false) => '0',
                        Some(true) => '1',
                    })
                    .collect();
                for (average, word) in averages.iter_mut().zip(OBSERVABLES) {
                    *average += basis_value(word, &output) / SHOTS as f64;
                }
            }
            for (average, word) in averages.into_iter().zip(OBSERVABLES) {
                let mut sum = LossyPauliSum::<f64>::new(2);
                sum += (word, 1.0);
                backward(&mut sum);
                let expected: f64 = sum
                    .iter()
                    .map(|(key, c)| c * basis_value(&key.to_string(), initial))
                    .sum();
                assert_abs_diff_eq!(average, expected, epsilon = 0.035);
            }
        }
    }
}

#[test]
fn pauli_probability_order_agrees_with_tableau_for_every_error() {
    for error in 0..15 {
        let mut p = [0.0; 15];
        p[error] = 1.0;
        for (word, plus) in [("ZI", false), ("IZ", false), ("XI", true), ("IX", true)] {
            let mut sim = Sim::new_with_seed(2, 0);
            if plus {
                sim.h(0);
                sim.h(1);
            }
            sim.two_qubit_pauli_error(0, 1, p);
            if plus {
                sim.h(0);
                sim.h(1);
            }
            let q = usize::from(word.starts_with('I'));
            let measured = if sim.measure(q).unwrap() { -1.0 } else { 1.0 };
            let mut sum = LossyPauliSum::<f64>::new(2);
            sum += (word, 1.0);
            sum.two_qubit_pauli_error(0, 1, p);
            assert_eq!(sum.iter().next().unwrap().1, measured);
        }
    }
}

/// Includes lost targets and lost partners: the tableau skips errors on lost
/// qubits but still applies the partner's half of a two-qubit Pauli.
#[test]
fn pauli_and_depolarizing_averages_agree_with_tableau() {
    let channels: [ChannelPair; 4] = [
        (
            |s| s.pauli_error(0, [0.1, 0.2, 0.3]),
            |s| s.pauli_error(0, [0.1, 0.2, 0.3]),
        ),
        (|s| s.depolarize1(0, 0.4), |s| s.depolarize1(0, 0.4)),
        (|s| s.depolarize2(0, 1, 0.4), |s| s.depolarize2(0, 1, 0.4)),
        (
            |s| s.two_qubit_pauli_error(0, 1, distinct_rates()),
            |s| s.two_qubit_pauli_error(0, 1, distinct_rates()),
        ),
    ];
    assert_matches_tableau(&["00", "11", "L1", "1L", "LL"], &channels);
}
