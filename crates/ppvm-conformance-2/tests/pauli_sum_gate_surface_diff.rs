// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness for the gate surface restored in the gate-parity
//! pass — `CliffordExtensions` (`s_dag`, `sqrt_x`, `sqrt_x_dag`, `sqrt_y`,
//! `sqrt_y_dag`, `cy`), `Projection` (`p0`/`p1`), `TwoQubitPauliError`,
//! `Depolarizing`, `Depolarizing2`, `AmplitudeDamping` and `RotXY::r` — against
//! the old `ppvm-pauli-sum::PauliSum<f64>` reference.
//!
//! Most tests build matched old/new sums, apply the identical gate on both, and
//! require term-for-term agreement. `Projection` is the documented exception:
//! the Lean oracle proves old's quadratic halving and `X`/`Y` arm wrong, so those
//! tests pin the deliberate Lean-correct divergence instead.

use ppvm_conformance_2::{
    assert_supports_match, build_new_sum, build_old_sum, new_support, old_support, random_terms,
    seeded_rng,
};

use ppvm_traits::traits::{
    AmplitudeDamping as OldAmplitudeDamping, CliffordExtensions as OldCliffordExtensions,
    Depolarizing as OldDepolarizing, Depolarizing2 as OldDepolarizing2,
    Projection as OldProjection, RotXY as OldRotXY, TwoQubitPauliError as OldTwoQubitPauliError,
};
use ppvm_traits_2::{
    AmplitudeDamping as NewAmplitudeDamping, CliffordExtensions as NewCliffordExtensions,
    Depolarizing as NewDepolarizing, Depolarizing2 as NewDepolarizing2,
    Projection as NewProjection, RotXY as NewRotXY, TwoQubitPauliError as NewTwoQubitPauliError,
};

/// Comparison tolerance: these are one- or two-gate applications, so the drift is
/// a few ulp.
const TOL: f64 = 1e-12;

const N: usize = 4;

fn matched(seed: u64, count: usize) -> (ppvm_conformance_2::OldSum, ppvm_conformance_2::NewSum) {
    let mut rng = seeded_rng(seed);
    let terms = random_terms(&mut rng, N, count);
    (build_old_sum(N, &terms), build_new_sum(N, &terms))
}

// --- CliffordExtensions ------------------------------------------------------

#[test]
fn clifford_extensions_match_old() {
    type OldGate = fn(&mut ppvm_conformance_2::OldSum, usize);
    type NewGate = fn(&mut ppvm_conformance_2::NewSum, usize);
    let gates: [(&str, OldGate, NewGate); 5] = [
        ("s_dag", |s, q| s.s_dag(q), |s, q| s.s_dag(q)),
        ("sqrt_x", |s, q| s.sqrt_x(q), |s, q| s.sqrt_x(q)),
        ("sqrt_x_dag", |s, q| s.sqrt_x_dag(q), |s, q| s.sqrt_x_dag(q)),
        ("sqrt_y", |s, q| s.sqrt_y(q), |s, q| s.sqrt_y(q)),
        ("sqrt_y_dag", |s, q| s.sqrt_y_dag(q), |s, q| s.sqrt_y_dag(q)),
    ];

    for (i, (name, old_gate, new_gate)) in gates.into_iter().enumerate() {
        for q in 0..N {
            let (mut old, mut new) = matched(0xC1FF + i as u64, 24);
            old_gate(&mut old, q);
            new_gate(&mut new, q);
            assert_supports_match(&old, &new, TOL);
            let _ = name;
        }
    }
}

#[test]
fn cy_matches_old() {
    for (a, b) in [(0usize, 1usize), (1, 3), (3, 0)] {
        let (mut old, mut new) = matched(0xC1FF_0002, 24);
        old.cy(a, b);
        new.cy(a, b);
        assert_supports_match(&old, &new, TOL);
    }
}

/// The single-qubit Cliffords are bit-level specializations in both crates; a
/// deep alternating sequence would catch a sign convention that only shows up
/// under composition.
#[test]
fn clifford_extension_sequence_matches_old() {
    let (mut old, mut new) = matched(0xC1FF_0003, 32);
    for q in 0..N {
        old.s_dag(q);
        new.s_dag(q);
        old.sqrt_x(q);
        new.sqrt_x(q);
        old.sqrt_y_dag(q);
        new.sqrt_y_dag(q);
        old.cy(q, (q + 1) % N);
        new.cy(q, (q + 1) % N);
        assert_supports_match(&old, &new, TOL);
    }
}

// --- RotXY -------------------------------------------------------------------

/// Contract 10: the sub-rotations are emitted in Heisenberg (backward) order, so
/// `r(q, π/2, θ) == ry(q, θ)`. Diffed against old, which pins the same identity.
#[test]
fn rot_xy_matches_old() {
    for &axis in &[0.0_f64, std::f64::consts::FRAC_PI_2, 1.3, -0.7] {
        for &theta in &[0.2_f64, 2.1] {
            let (mut old, mut new) = matched(0x8000, 16);
            old.r(0, axis, theta);
            new.r(0, axis, theta);
            assert_supports_match(&old, &new, TOL);
        }
    }
}

// --- Diagonal channels -------------------------------------------------------

#[test]
fn two_qubit_pauli_error_matches_old() {
    // A *mixed* probability vector, not a one-hot: a transposed index in the
    // hand-written anticommuting-set tables is invisible on one-hot inputs.
    let p: [f64; 15] = std::array::from_fn(|i| 0.001 * (i as f64 + 1.0));
    for (a, b) in [(0usize, 1usize), (2, 3), (0, 3)] {
        let (mut old, mut new) = matched(0x2100_u64.wrapping_add(a as u64), 40);
        old.two_qubit_pauli_error(a, b, p);
        new.two_qubit_pauli_error(a, b, p);
        assert_supports_match(&old, &new, TOL);
    }
}

#[test]
fn depolarizing_matches_old() {
    for &prob in &[0.0_f64, 0.05, 0.75] {
        for q in 0..N {
            let (mut old, mut new) = matched(0xDE00, 32);
            old.depolarize1(q, prob);
            new.depolarize1(q, prob);
            assert_supports_match(&old, &new, TOL);
        }
    }
}

#[test]
fn depolarizing2_matches_old() {
    for &prob in &[0.0_f64, 0.05, 0.9375] {
        for (a, b) in [(0usize, 1usize), (1, 2), (0, 3)] {
            let (mut old, mut new) = matched(0xDE02, 32);
            old.depolarize2(a, b, prob);
            new.depolarize2(a, b, prob);
            assert_supports_match(&old, &new, TOL);
        }
    }
}

/// Amplitude damping is the one branching channel here; its `Z → I` branch must
/// accumulate onto an existing `I` on both sides (contract 3(c)).
#[test]
fn amplitude_damping_matches_old() {
    for &gamma in &[0.0_f64, 0.25, 1.0] {
        for q in 0..N {
            let (mut old, mut new) = matched(0xADAA, 40);
            old.amplitude_damping(q, gamma);
            new.amplitude_damping(q, gamma);
            assert_supports_match(&old, &new, TOL);
        }
    }
}

// --- Projection --------------------------------------------------------------

/// Old computes `c²/2`; Lean proves the projector is linear and therefore uses
/// `c/2`. Pin both the divergence and the corrected value.
#[test]
fn projection_uses_the_lean_correct_linear_halving() {
    let terms = vec![("I".to_string(), 2.0)];
    let mut old = build_old_sum(1, &terms);
    let mut new = build_new_sum(1, &terms);
    old.p0(0);
    new.p0(0);

    assert_eq!(
        old_support(&old),
        vec![("I".into(), 2.0), ("Z".into(), 2.0)]
    );
    assert_eq!(
        new_support(&new),
        vec![("I".into(), 1.0), ("Z".into(), 1.0)]
    );

    // A genuine projector is idempotent.
    new.p0(0);
    assert_eq!(
        new_support(&new),
        vec![("I".into(), 1.0), ("Z".into(), 1.0)]
    );
}

/// The matrix oracle gives `ΠXΠ = ΠYΠ = 0`; old leaves both untouched. Zero
/// entries remain present until explicit reduction.
#[test]
fn projection_annihilates_x_and_y_as_the_matrix_oracle_requires() {
    let terms = vec![("X".to_string(), 2.0), ("Y".to_string(), -3.0)];
    let mut old = build_old_sum(1, &terms);
    let mut new = build_new_sum(1, &terms);
    old.p0(0);
    new.p0(0);

    assert_eq!(old_support(&old), terms);
    assert_eq!(
        new_support(&new),
        vec![("X".into(), 0.0), ("Y".into(), 0.0)]
    );
    new.reduce();
    assert!(new.is_empty());
}
