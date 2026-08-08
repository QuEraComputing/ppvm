// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Behaviour tests for the gate surface restored in this pass: the native
//! two-qubit rotations (`RotationTwo`), `RotXY`, `CliffordExtensions`,
//! `Projection`, the remaining noise channels, and the `Display`/`Debug`
//! rendering.
//!
//! The parity tests are ported from the old crate's own: `rot2.rs`'s
//! `rxx_matches_generic` / `ryy_matches_generic` / `rzz_matches_generic` and
//! `rzz_explicit_values` / `rzz_non_adjacent_addressing`, and `rot1.rs`'s
//! `test_r`.

use ppvm_pauli_sum_2::{PauliSum, PauliWord};
use ppvm_traits_2::{
    AmplitudeDamping, Clifford, CliffordExtensions, Depolarizing, Depolarizing2, PauliBits,
    Projection, RotXY, RotationOne, RotationTwo, TwoQubitPauliError, Word,
};

fn pw(s: &str) -> PauliWord {
    PauliWord::from(s)
}

fn sum_with(n: usize, word: &str) -> PauliSum {
    PauliSum::from_terms(n, [(pw(word), 1.0)])
}

fn approx(sum: &PauliSum, key: &str, want: f64) {
    let got = sum.get(&pw(key)).unwrap_or(0.0);
    assert!(
        (got - want).abs() < 1e-12,
        "coeff at {key}: got {got}, want {want}"
    );
}

/// Exact-support equality: same key set, same coefficients (to 1e-12).
fn assert_same(got: &PauliSum, want: &PauliSum, ctx: &str) {
    assert_eq!(
        got.len(),
        want.len(),
        "{ctx}: len differs\n{got:?}\n{want:?}"
    );
    for (k, v) in got.iter() {
        let w = want.get(&k).unwrap_or_else(|| panic!("{ctx}: {k} missing"));
        assert!((v - w).abs() < 1e-12, "{ctx}: {k}: got {v}, want {w}");
    }
}

const PAULI_CHARS: [char; 4] = ['I', 'X', 'Y', 'Z'];

// --- RotationTwo: each fast path against the generic `rotate_2`. --------------

fn assert_matches_generic(axis: [u8; 2], special: impl Fn(&mut PauliSum, usize, usize, f64)) {
    let theta = 0.7_f64;
    for &p in &PAULI_CHARS {
        for &q in &PAULI_CHARS {
            let word = format!("{p}{q}");

            let mut got = sum_with(2, &word);
            special(&mut got, 0, 1, theta);

            let mut want = sum_with(2, &word);
            want.rotate_2(axis, axis, 0, 1, theta);

            assert_same(&got, &want, &format!("fast path vs rotate_2 for {word}"));
        }
    }
}

#[test]
fn rxx_matches_generic() {
    assert_matches_generic([1, 0], |s, a, b, t| s.rxx(a, b, t));
}

#[test]
fn ryy_matches_generic() {
    assert_matches_generic([1, 1], |s, a, b, t| s.ryy(a, b, t));
}

#[test]
fn rzz_matches_generic() {
    assert_matches_generic([0, 1], |s, a, b, t| s.rzz(a, b, t));
}

/// Hand-computed values independent of `rotate_2`/`comm_2`, so a bug shared by
/// both paths cannot hide (old's `rzz_explicit_values`).
#[test]
fn rzz_explicit_values() {
    let t = 0.7_f64;

    // X₀I₁ --rzz--> cos·XI − sin·YZ (anticommutes; the X carrier gives −1).
    let mut s = sum_with(2, "XI");
    s.rzz(0, 1, t);
    approx(&s, "XI", t.cos());
    approx(&s, "YZ", -t.sin());
    assert_eq!(s.len(), 2);

    // Y₀I₁ --rzz--> cos·YI + sin·XZ (the Y carrier gives +1).
    let mut s = sum_with(2, "YI");
    s.rzz(0, 1, t);
    approx(&s, "YI", t.cos());
    approx(&s, "XZ", t.sin());

    // ZZ commutes with the ZZ generator → unchanged.
    let mut s = sum_with(2, "ZZ");
    s.rzz(0, 1, t);
    assert_eq!(s.len(), 1);
    approx(&s, "ZZ", 1.0);
}

/// A diagonal rotation on non-adjacent qubits must address the right slots
/// (old's `rzz_non_adjacent_addressing`).
#[test]
fn rzz_non_adjacent_addressing() {
    let theta = 0.7_f64;
    for &p in &PAULI_CHARS {
        for &q in &PAULI_CHARS {
            let word = format!("{p}I{q}");

            let mut got = sum_with(3, &word);
            got.rzz(0, 2, theta);

            let mut want = sum_with(3, &word);
            want.rotate_2([0, 1], [0, 1], 0, 2, theta);

            assert_same(&got, &want, &format!("rzz(0,2) vs rotate_2 for {word}"));
        }
    }
}

/// The native `rzz` must agree with the `cnot; rz; cnot` decomposition it
/// replaces — the identity the headline integration bench relied on before the
/// native kernel existed.
#[test]
fn rzz_matches_cnot_rz_cnot_decomposition() {
    let theta = 0.37_f64;
    for &p in &PAULI_CHARS {
        for &q in &PAULI_CHARS {
            let word = format!("{p}{q}");

            let mut native = sum_with(2, &word);
            native.rzz(0, 1, theta);

            let mut decomposed = sum_with(2, &word);
            decomposed.cnot(0, 1);
            decomposed.rz(1, theta);
            decomposed.cnot(0, 1);

            assert_same(
                &native,
                &decomposed,
                &format!("native rzz vs cnot;rz;cnot for {word}"),
            );
        }
    }
}

/// Non-bit axis components are rejected. Old's guard is `> 3` with a
/// "cannot be L" message, which lets `2`/`3` through into an out-of-bounds table
/// index (suspected old bug 5); the `-2` bound is `> 1`.
#[test]
#[should_panic(expected = "rotation axis components must be 0 or 1")]
fn rotate_2_rejects_non_bit_axis() {
    let mut s = sum_with(2, "XX");
    s.rotate_2([2, 0], [0, 1], 0, 1, 0.5);
}

/// `rotate_2` must not truncate or reduce: the identity angle still merges the
/// exactly-zero branch (old's `map_insert` `add_assign`s every produced term).
#[test]
fn rzz_zero_angle_still_inserts_the_zero_branch() {
    let mut s = sum_with(2, "XI");
    s.rzz(0, 1, 0.0);
    assert_eq!(s.len(), 2, "the sin = 0 branch must still be inserted");
    approx(&s, "XI", 1.0);
    approx(&s, "YZ", 0.0);
}

// --- RotXY: the Heisenberg (backward) sub-rotation order. --------------------

/// `r(q, π/2, θ) == ry(q, θ)` and **not** `ry(q, −θ)`: the case that
/// distinguishes the Heisenberg order from the Schrödinger one (old's
/// `rot1.rs::test_r`).
#[test]
fn r_is_heisenberg_ordered() {
    use std::f64::consts::FRAC_PI_2;
    let theta = 2.1_f64;

    let mut via_r = sum_with(1, "Z");
    via_r.r(0, 0.0, theta);
    let mut via_rx = sum_with(1, "Z");
    via_rx.rx(0, theta);
    assert!((via_r.overlap(&via_rx) - 1.0).abs() < 1e-9);

    let mut via_r = sum_with(1, "Z");
    via_r.r(0, FRAC_PI_2, theta);
    let mut via_ry = sum_with(1, "Z");
    via_ry.ry(0, theta);
    assert!(
        (via_r.overlap(&via_ry) - 1.0).abs() < 1e-9,
        "r(π/2, θ) must equal ry(θ), not ry(−θ)"
    );
}

// --- CliffordExtensions: the conjugation tables. -----------------------------

/// Conjugate the one-qubit Paulis and read back `(sign, image)`.
fn conj1(gate: impl Fn(&mut PauliSum, usize), input: &str) -> (f64, String) {
    let mut s = sum_with(1, input);
    gate(&mut s, 0);
    assert_eq!(s.len(), 1, "a single-qubit Clifford is a bijection");
    let (k, v) = s.iter().next().unwrap();
    (v, k.to_string())
}

/// One row of the single-qubit conjugation table: the gate's name, the gate
/// itself, and its `(sign, image)` action on `X`, `Y`, `Z` in that order.
type Conj1Case = (
    &'static str,
    fn(&mut PauliSum, usize),
    [(f64, &'static str); 3],
);

#[test]
fn clifford_extension_tables_match_the_trait_doc() {
    // (gate, [X image, Y image, Z image]) — the table in `CliffordExtensions`'
    // docs, which is old's `sum/clifford.rs` behaviour.
    let cases: [Conj1Case; 6] = [
        ("s", |s, q| s.s(q), [(-1.0, "Y"), (1.0, "X"), (1.0, "Z")]),
        (
            "s_dag",
            |s, q| s.s_dag(q),
            [(1.0, "Y"), (-1.0, "X"), (1.0, "Z")],
        ),
        (
            "sqrt_x",
            |s, q| s.sqrt_x(q),
            [(1.0, "X"), (-1.0, "Z"), (1.0, "Y")],
        ),
        (
            "sqrt_x_dag",
            |s, q| s.sqrt_x_dag(q),
            [(1.0, "X"), (1.0, "Z"), (-1.0, "Y")],
        ),
        (
            "sqrt_y",
            |s, q| s.sqrt_y(q),
            [(1.0, "Z"), (1.0, "Y"), (-1.0, "X")],
        ),
        (
            "sqrt_y_dag",
            |s, q| s.sqrt_y_dag(q),
            [(-1.0, "Z"), (1.0, "Y"), (1.0, "X")],
        ),
    ];

    for (name, gate, table) in cases {
        for (input, (want_sign, want_image)) in ["X", "Y", "Z"].iter().zip(table) {
            let (sign, image) = conj1(gate, input);
            assert_eq!(image, want_image, "{name} on {input}: image");
            assert_eq!(sign, want_sign, "{name} on {input}: sign");
        }
        // The identity is fixed by every Clifford.
        assert_eq!(conj1(gate, "I"), (1.0, "I".to_string()), "{name} on I");
    }
}

/// The bit-level specializations of `h`/`s` must agree with the phased-word
/// round-trip they replaced: `h` is its own inverse and `s·s = z` up to sign.
#[test]
fn h_and_s_specializations_are_consistent() {
    for &p in &PAULI_CHARS {
        let word = p.to_string();

        // H is an involution on the sum.
        let mut s = sum_with(1, &word);
        s.h(0);
        s.h(0);
        assert_same(&s, &sum_with(1, &word), &format!("h∘h on {word}"));

        // S∘S† is the identity.
        let mut s = sum_with(1, &word);
        s.s(0);
        s.s_dag(0);
        assert_same(&s, &sum_with(1, &word), &format!("s∘s† on {word}"));

        // √X∘√X = X conjugation (a pure sign).
        let mut got = sum_with(1, &word);
        got.sqrt_x(0);
        got.sqrt_x(0);
        let mut want = sum_with(1, &word);
        want.x(0);
        assert_same(&got, &want, &format!("√X∘√X vs x on {word}"));
    }
}

/// `cy` still goes through the phased word (old's macro); check its table.
#[test]
fn cy_conjugation_table() {
    // CY: X₀ ↦ X₀Y₁.
    let mut s = sum_with(2, "XI");
    s.cy(0, 1);
    assert_eq!(s.len(), 1);
    approx(&s, "XY", 1.0);
    // Z₀ is fixed by a control-side Z.
    let mut s = sum_with(2, "ZI");
    s.cy(0, 1);
    approx(&s, "ZI", 1.0);
}

// --- Diagonal noise: the key set is invariant (contract 4). ------------------

#[test]
fn diagonal_channels_never_insert_or_remove() {
    let seed: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("II"), 1.0),
            (pw("XI"), 1.0),
            (pw("YZ"), 1.0),
            (pw("ZY"), 1.0),
        ],
    );

    let mut s = seed.clone();
    s.two_qubit_pauli_error(0, 1, [0.01; 15]);
    assert_eq!(s.len(), seed.len());
    for (k, _) in seed.iter() {
        assert!(s.contains_key(&k), "{k} vanished");
    }

    let mut s = seed.clone();
    s.depolarize1(0, 0.75); // factor 1 − 4·0.75/3 == 0.0 exactly
    assert_eq!(s.len(), seed.len(), "a zero factor must not remove terms");
    approx(&s, "XI", 0.0);
    approx(&s, "II", 1.0);

    let mut s = seed.clone();
    s.depolarize2(0, 1, 0.1);
    assert_eq!(s.len(), seed.len());
    let f = 1.0 - 0.1 * (16.0 / 15.0);
    approx(&s, "II", 1.0);
    approx(&s, "XI", f);
    approx(&s, "YZ", f);
}

/// A one-hot two-qubit probability vector reproduces the single-qubit channel it
/// degenerates to. `p[ZI] = p[11]` flips every term that anticommutes with `Z₀`.
#[test]
fn two_qubit_pauli_error_one_hot_zi() {
    let mut p = [0.0_f64; 15];
    p[11] = 0.1; // ZI
    let mut s: PauliSum = PauliSum::from_terms(
        2,
        [
            (pw("II"), 1.0),
            (pw("XI"), 1.0),
            (pw("ZI"), 1.0),
            (pw("IX"), 1.0),
        ],
    );
    s.two_qubit_pauli_error(0, 1, p);
    approx(&s, "II", 1.0); // commutes
    approx(&s, "XI", 1.0 - 0.2); // anticommutes with Z₀
    approx(&s, "ZI", 1.0); // commutes
    approx(&s, "IX", 1.0); // commutes with Z₀
}

/// Amplitude damping's `Z → I` branch must **accumulate** onto an existing `I`,
/// never overwrite it (behavioural contract 3(c)).
#[test]
fn amplitude_damping_branch_accumulates_onto_identity() {
    let gamma = 0.25_f64;
    let mut s: PauliSum = PauliSum::from_terms(1, [(pw("I"), 2.0), (pw("Z"), 4.0)]);
    s.amplitude_damping(0, gamma);
    assert_eq!(s.len(), 2);
    approx(&s, "I", 2.0 + gamma * 4.0);
    approx(&s, "Z", 4.0 * (1.0 - gamma));

    // X / Y damp by √(1−γ) with no branch.
    let mut s: PauliSum = PauliSum::from_terms(1, [(pw("X"), 1.0), (pw("Y"), 1.0)]);
    s.amplitude_damping(0, gamma);
    assert_eq!(s.len(), 2);
    approx(&s, "X", (1.0 - gamma).sqrt());
    approx(&s, "Y", (1.0 - gamma).sqrt());
}

// --- Projection -------------------------------------------------------------

/// `p0`/`p1` use the Lean-correct linear projector, including outside old's
/// unit-coefficient blind spot.
#[test]
fn projection_splits_identity_and_z() {
    let mut s = sum_with(1, "I");
    s.p0(0);
    assert_eq!(s.len(), 2);
    approx(&s, "I", 0.5);
    approx(&s, "Z", 0.5);

    let mut s = sum_with(1, "I");
    s.p1(0);
    approx(&s, "I", 0.5);
    approx(&s, "Z", -0.5);

    // A non-unit coefficient distinguishes the linear `c/2` map from old's
    // quadratic `c²/2`.
    let mut s: PauliSum = PauliSum::from_terms(1, [(pw("I"), 2.0)]);
    s.p0(0);
    approx(&s, "I", 1.0);
    approx(&s, "Z", 1.0);
    s.p0(0);
    approx(&s, "I", 1.0);
    approx(&s, "Z", 1.0);

    // X / Y are annihilated. The zero key remains until explicit `reduce()`,
    // matching the engine's no-implicit-reduce contract.
    let mut s = sum_with(1, "X");
    s.p0(0);
    assert_eq!(s.len(), 1);
    approx(&s, "X", 0.0);
    s.reduce();
    assert!(s.is_empty());
}

// --- Display / Debug (contract 13) ------------------------------------------

#[test]
fn display_and_debug_render_weight_sorted_terms() {
    let s: PauliSum = PauliSum::from_terms(2, [(pw("XY"), 0.5)]);
    assert_eq!(format!("{s}"), "0.500 * XY");
    assert_eq!(format!("{s:?}"), "0.50000000 * XY");

    // Ascending Pauli weight, joined by " + ".
    let s: PauliSum = PauliSum::from_terms(2, [(pw("XY"), 1.0), (pw("IZ"), 2.0)]);
    let rendered = format!("{s}");
    assert_eq!(rendered, "2.000 * IZ + 1.000 * XY", "weight order");

    // An empty sum renders as the empty string.
    let s: PauliSum = PauliSum::new(2);
    assert_eq!(format!("{s}"), "");
}

/// `weight()` really is the sort key (not the map order).
#[test]
fn display_sorts_by_pauli_weight() {
    let s: PauliSum = PauliSum::from_terms(3, [(pw("XXX"), 1.0), (pw("IIZ"), 1.0)]);
    assert_eq!(pw("IIZ").weight(), 1);
    assert_eq!(pw("XXX").weight(), 3);
    assert!(format!("{s}").starts_with("1.000 * IIZ"));
}

// --- The non-lossy word keeps every loss branch dead (feature 11). -----------

#[test]
fn plain_pauli_word_is_never_lost() {
    let w = pw("XYZ");
    for q in 0..3 {
        assert!(!w.is_lost(q));
    }
    assert_eq!(w.loss_weight(), 0);
}
