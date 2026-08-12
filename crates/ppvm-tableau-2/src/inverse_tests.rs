// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Tests for the inverse-tableau signs ([`crate::inverse`]).
//!
//! Two independent oracles, because the sign rules are a derivation rather than
//! a transcription:
//!
//! 1. `U·(U†PU)·U† = P` — [`Tableau::assert_inverse_consistent`], which checks
//!    every rule against the *definition* of the inverse rather than against a
//!    rearrangement of the formula it feeds.
//! 2. The decomposition phase and the deterministic outcome against the row
//!    folds they replace, which are the shipped, differentially-tested paths.

use ppvm_traits_2::{Clifford, CliffordExtensions, Measure, Pauli};
use rand::{RngExt, SeedableRng, rngs::SmallRng};

use crate::U256;
use crate::data::{GeneralizedTableau, Tableau};

/// Every gate the crate has, applied to a fresh frame and to an already
/// scrambled one — one gate per case, so a failure names the rule.
#[test]
fn every_gate_rule_holds_in_isolation() {
    for scramble in [false, true] {
        for gate in 0..GATE_COUNT {
            let mut tab = Tableau::<fxhash::FxBuildHasher>::new(5);
            if scramble {
                scramble_frame(&mut tab, &mut SmallRng::seed_from_u64(7), 12);
            }
            apply_indexed_gate(&mut tab, gate, 1, 3);
            tab.assert_inverse_consistent();
        }
    }
}

/// Random Clifford circuits, checked after **every** gate, at widths that
/// straddle the word and block boundaries.
#[test]
fn random_circuits_keep_the_inverse_signs() {
    for (n, seed) in [(1usize, 0u64), (2, 1), (7, 2), (64, 3), (65, 4), (100, 5)] {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut tab = Tableau::<fxhash::FxBuildHasher>::new(n);
        for _ in 0..40 {
            let (a, b) = distinct_pair(n, &mut rng);
            let limit = if n == 1 { GATES_1 } else { GATE_COUNT };
            apply_indexed_gate(&mut tab, rng.random_range(0..limit), a, b);
            tab.assert_inverse_consistent();
        }
    }
}

/// The `O(1)` decomposition phase equals the fold, at every site and for every
/// Pauli — including the `Y` site rule, which no measurement path reaches.
#[test]
fn decomposition_phase_matches_the_fold() {
    for (n, seed) in [(1usize, 10u64), (3, 11), (64, 12), (70, 13)] {
        let mut tab = GeneralizedTableau::<U256>::new(n, 1e-12);
        scramble_frame(&mut tab.tableau, &mut SmallRng::seed_from_u64(seed), 60);
        for q in 0..n {
            for pauli in [Pauli::X, Pauli::Y, Pauli::Z] {
                let inverse = tab.compute_decomposition(q, pauli);
                tab.tableau.invalidate_inverse();
                let fold = tab.compute_decomposition(q, pauli);
                tab.tableau.rebuild_inverse_signs();
                assert_eq!(inverse, fold, "n={n} q={q} {pauli:?}");
            }
        }
    }
}

/// A deterministic outcome is one inverse row's sign. Checked against the
/// stabilizer fold on a frame that keeps every qubit in the Z basis, so every
/// measurement stays in case b.
#[test]
fn deterministic_outcome_matches_the_fold() {
    let n = 8;
    let mut rng = SmallRng::seed_from_u64(21);
    let mut tab = Tableau::<fxhash::FxBuildHasher>::new(n);
    for _ in 0..30 {
        let (a, b) = distinct_pair(n, &mut rng);
        match rng.random_range(0..3) {
            0 => Clifford::x(&mut tab, a),
            1 => Clifford::cnot(&mut tab, a, b),
            _ => Clifford::cz(&mut tab, a, b),
        }
        for q in 0..n {
            assert!(
                tab.find_z_anticommuting_stabilizer(q).is_none(),
                "the fixture must stay in case b"
            );
            let inverse = tab.get_deterministic_outcome(q);
            tab.invalidate_inverse();
            let fold = tab.get_deterministic_outcome(q);
            tab.rebuild_inverse_signs();
            assert_eq!(inverse, fold, "q={q}");
        }
    }
}

/// A case-a projection carries the signs through, on frames scrambled between
/// measurements so the pivot's neighbourhood is not the same twice.
///
/// The scramble is what makes the pivot columns *dense*, and density is what
/// [`blocks::prefer_gather`](crate::storage::blocks) branches on: at `n = 64` the
/// projection takes the transpose 60 times and gathers 93 (counted by
/// instrumenting the predicate), while at `n = 1`, `6` and `70` it always
/// gathers. Both branches are therefore the projection's own choice rather than
/// a caller's guard. A layer of `H` does *not* do this — it leaves one X per
/// column, which is the sparse side of the predicate at every width.
#[test]
fn projection_keeps_the_inverse_signs() {
    for (n, seed) in [(1usize, 30u64), (6, 31), (64, 32), (70, 33)] {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mut tab = GeneralizedTableau::<U256>::new(n, 1e-12);
        for round in 0..3 {
            scramble_frame(&mut tab.tableau, &mut rng, 8 * n);
            for q in 0..n {
                tab.measure(q, &mut rng);
                tab.tableau.assert_inverse_consistent();
            }
            assert!(
                tab.tableau.inverse_valid(),
                "n={n} round={round}: the signs must survive a projection"
            );
        }
    }
}

/// A projection under an outer row guard — `measure_all`'s shape, where the
/// inverse update's site reads are `memcpy`s and the pivot's own selectors are
/// strided.
#[test]
fn projection_keeps_the_signs_under_a_row_guard() {
    let n = 12;
    let mut rng = SmallRng::seed_from_u64(41);
    let mut tab = GeneralizedTableau::<U256>::new(n, 1e-12);
    for q in 0..n {
        Clifford::h(&mut tab, q);
    }
    for q in 0..n.saturating_sub(1) {
        Clifford::cnot(&mut tab, q, q + 1);
    }
    let record = tab.measure_all(&mut rng);
    assert_eq!(record.len(), n);
    tab.tableau.assert_inverse_consistent();
}

/// The primitives that move bits and signs independently are not Cliffords, so
/// they must abandon the signs rather than leave them stale.
#[test]
fn non_clifford_primitives_abandon_the_signs() {
    use ppvm_traits_2::{PhaseTrack, StabilizerFrame, SymplecticColumns};

    /// A named frame mutation with no inverse rule.
    type Case = (&'static str, fn(&mut Tableau));

    let cases: [Case; 4] = [
        ("swap_xz", |t| t.swap_xz(1)),
        ("xor_z_from_x", |t| t.xor_z_from_x(1)),
        ("flip_phase_where_xz", |t| t.flip_phase_where_xz(1)),
        ("row_multiply", |t| t.row_multiply(0, 3)),
    ];
    for (name, op) in cases {
        let mut tab = Tableau::<fxhash::FxBuildHasher>::new(4);
        assert!(tab.inverse_valid());
        op(&mut tab);
        assert!(!tab.inverse_valid(), "{name} left the signs in place");
    }
}

// ─── Helpers ──────────────────────────────────────────────────────────────

/// How many one-qubit gates [`apply_indexed_gate`] knows; the rest are pairs.
const GATES_1: usize = 10;
/// Every gate it knows.
const GATE_COUNT: usize = GATES_1 + 3;

/// Apply gate `i` of the crate's whole Clifford set.
fn apply_indexed_gate<H>(tab: &mut Tableau<H>, i: usize, a: usize, b: usize) {
    match i {
        0 => Clifford::x(tab, a),
        1 => Clifford::y(tab, a),
        2 => Clifford::z(tab, a),
        3 => Clifford::h(tab, a),
        4 => Clifford::s(tab, a),
        5 => CliffordExtensions::s_dag(tab, a),
        6 => CliffordExtensions::sqrt_x(tab, a),
        7 => CliffordExtensions::sqrt_x_dag(tab, a),
        8 => CliffordExtensions::sqrt_y(tab, a),
        9 => CliffordExtensions::sqrt_y_dag(tab, a),
        10 => Clifford::cnot(tab, a, b),
        11 => Clifford::cz(tab, a, b),
        _ => CliffordExtensions::cy(tab, a, b),
    }
}

/// Two distinct qubits, or `(0, 0)` on a one-qubit frame — where the caller must
/// pick a one-qubit gate.
fn distinct_pair(n: usize, rng: &mut SmallRng) -> (usize, usize) {
    if n == 1 {
        return (0, 0);
    }
    let a = rng.random_range(0..n);
    let mut b = rng.random_range(0..n);
    while b == a {
        b = rng.random_range(0..n);
    }
    (a, b)
}

/// A pseudo-random Clifford circuit.
fn scramble_frame<H>(tab: &mut Tableau<H>, rng: &mut SmallRng, gates: usize) {
    let n = tab.n_qubits();
    for _ in 0..gates {
        let (a, b) = distinct_pair(n, rng);
        let limit = if n == 1 { GATES_1 } else { GATE_COUNT };
        apply_indexed_gate(tab, rng.random_range(0..limit), a, b);
    }
}
