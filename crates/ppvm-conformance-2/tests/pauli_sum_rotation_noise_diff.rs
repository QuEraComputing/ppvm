// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential correctness for the **non-Clifford** surface of
//! `ppvm-pauli-sum-2::PauliSum<f64>` — the branching single-qubit rotations
//! (`rx`/`ry`/`rz`) and the diagonal `pauli_error` channel — against the old
//! `ppvm-pauli-sum::PauliSum<f64>` reference, driven by the shared seeded
//! generators of `ppvm-conformance-2` and reusing its `build_old_sum`/
//! `build_new_sum` harness pair.
//!
//! What each test exercises (design: `traits-2-configuration-and-hashing.md`
//! §"Behavioral traits", Lean `Instantiations/Rotation.lean` + `Algebra/Noise.lean`):
//!
//! * **(a) rotation replay** — a seeded random `rx`/`ry`/`rz` circuit (random
//!   angles, random qubits) on matched old and new sums; after *every* gate the
//!   full support + coefficients must agree. A rotation is a genuine fan-out
//!   (1→2 terms) whose `iGP` branch may collide with an existing key, so this
//!   drives the branch producer **and** the collision-merge in `accumulate_batch`.
//! * **(b) exhaustive single-qubit** — `rx`/`ry`/`rz` by a fixed angle on each of
//!   `I`/`X`/`Y`/`Z`, checked against an *independent* hand-derived `cos`/`sin`
//!   reference table (the same values `ppvm-pauli-sum`'s own unit tests assert),
//!   so a regression in *both* crates cannot hide behind the diff.
//! * **(c) `pauli_error`** — random probability triples on random qubits vs old,
//!   asserting the per-term `λ_P` scaling matches term for term.
//! * **(d) collision merge** — two rotation branches that land on the same key,
//!   checked to add correctly in both crates and against the closed form.
//!
//! We compare **observable algebra, not raw hash digests** (the crates' hash
//! finalization folds differ by design).

use ppvm_conformance_2::{
    NewSum, OldSum, assert_close, build_new_sum, build_old_sum, new_support, old_support,
    random_terms, seeded_rng,
};

// Rotation + noise traits, aliased to keep the old and new namespaces apart.
use ppvm_traits::traits::{PauliError as OldPauliError, RotationOne as OldRotationOne};
use ppvm_traits_2::{PauliError as NewPauliError, RotationOne as NewRotationOne};

use rand::RngExt;
use rand::rngs::StdRng;

/// Seeds every property sweep replays under.
const SEEDS: [u64; 10] = [1, 2, 3, 7, 42, 99, 123, 777, 2024, 31337];
/// Qubit widths exercised (≤ 64, the shared backing capacity).
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 8, 12];
/// Coefficient comparison tolerance for the `f64` backends. Loosened from the
/// Clifford suite's `1e-9`: a rotation replay accumulates `sin`/`cos` products
/// (and the two crates merge colliding keys in different iteration orders), so a
/// few ulp of drift per gate is expected.
const TOL: f64 = 1e-7;
/// Coefficients whose magnitude falls below this are treated as *absent*. A
/// rotation merge can drive a coefficient to a physical near-zero; the new crate
/// drops only the *exact* zeros its `reduce` sees, while the old crate leaves the
/// residue in its map, so we compare on the same above-floor support in both.
/// Well below any generated `O(1)` coefficient, well above `f64` merge noise.
const FLOOR: f64 = 1e-9;

/// Assert the old and new supports agree after filtering out below-`FLOOR`
/// residue in both: same above-floor keys in canonical order, coefficients close.
#[track_caller]
fn assert_rot_supports_match(old: &OldSum, new: &NewSum, tol: f64) {
    let filt = |v: Vec<(String, f64)>| -> Vec<(String, f64)> {
        v.into_iter().filter(|(_, c)| c.abs() > FLOOR).collect()
    };
    let os = filt(old_support(old));
    let ns = filt(new_support(new));
    assert_eq!(
        os.len(),
        ns.len(),
        "above-floor support size differs: old {} vs new {}\nold={os:?}\nnew={ns:?}",
        os.len(),
        ns.len()
    );
    for (o, n) in os.iter().zip(ns.iter()) {
        assert_eq!(o.0, n.0, "support key differs: old {} vs new {}", o.0, n.0);
        assert_close(o.1, n.1, tol.max(o.1.abs() * 1e-9));
    }
}

// ---------------------------------------------------------------------------
// (a) THE CORE — rotation replay with per-gate support+coefficient agreement.
// ---------------------------------------------------------------------------

/// A replayable single-qubit rotation gate — `rx`/`ry`/`rz` with an angle.
#[derive(Clone, Copy, Debug)]
enum Rot {
    Rx(usize, f64),
    Ry(usize, f64),
    Rz(usize, f64),
}

fn apply_old_rot(s: &mut OldSum, g: Rot) {
    match g {
        Rot::Rx(q, t) => s.rx(q, t),
        Rot::Ry(q, t) => s.ry(q, t),
        Rot::Rz(q, t) => s.rz(q, t),
    }
}

fn apply_new_rot(s: &mut NewSum, g: Rot) {
    match g {
        Rot::Rx(q, t) => s.rx(q, t),
        Rot::Ry(q, t) => s.ry(q, t),
        Rot::Rz(q, t) => s.rz(q, t),
    }
}

/// A random `rx`/`ry`/`rz` gate list on `n` qubits, `len` gates long, angles in
/// `(−π, π)`.
fn random_rot_ops(rng: &mut StdRng, n: usize, len: usize) -> Vec<Rot> {
    let angle = |rng: &mut StdRng| rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
    (0..len)
        .map(|_| {
            let q = rng.random_range(0..n);
            match rng.random_range(0..3usize) {
                0 => Rot::Rx(q, angle(rng)),
                1 => Rot::Ry(q, angle(rng)),
                _ => Rot::Rz(q, angle(rng)),
            }
        })
        .collect()
}

#[test]
fn rotation_replay_matches_old_per_gate() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            // A moderate initial support so branch fan-out and collision-merges
            // have something to act on. A rotation fans out (1→2 terms), so the
            // support grows toward the full `4ⁿ` basis; cap the gate count on wide
            // registers to keep the per-gate O(m log m) compare bounded while still
            // driving many branches and merges.
            let terms = random_terms(&mut rng, n, 12);
            let mut old = build_old_sum(n, &terms);
            let mut new = build_new_sum(n, &terms);

            // Agreement before any gate.
            assert_rot_supports_match(&old, &new, TOL);

            let gates = if n <= 5 { 30 } else { 12 };
            let ops = random_rot_ops(&mut rng, n, gates);
            for g in ops {
                apply_old_rot(&mut old, g);
                apply_new_rot(&mut new, g);
                assert_rot_supports_match(&old, &new, TOL);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// (b) EXHAUSTIVE single-qubit against an independent cos/sin reference.
// ---------------------------------------------------------------------------

/// Reference conjugation `g P g†` for a single-qubit rotation in the Heisenberg
/// picture, as a list of `(word, coeff)` terms. These are the *independent*
/// hand-derived values (identical to the ones `ppvm-pauli-sum`'s own `rot1`
/// unit tests assert): a genuine oracle, not a call into either crate.
///
///   rx: Y → cosθ·Y − sinθ·Z,  Z → cosθ·Z + sinθ·Y
///   ry: X → cosθ·X + sinθ·Z,  Z → cosθ·Z − sinθ·X
///   rz: X → cosθ·X − sinθ·Y,  Y → cosθ·Y + sinθ·X
fn reference_single(axis: char, p: &str, c: f64, theta: f64) -> Vec<(String, f64)> {
    let (s, k) = (theta.sin(), theta.cos());
    let two =
        |w0: &str, c0: f64, w1: &str, c1: f64| vec![(w0.to_string(), c0), (w1.to_string(), c1)];
    match (axis, p) {
        ('x', "I") | ('x', "X") => vec![(p.to_string(), c)],
        ('x', "Y") => two("Y", c * k, "Z", -c * s),
        ('x', "Z") => two("Z", c * k, "Y", c * s),
        ('y', "I") | ('y', "Y") => vec![(p.to_string(), c)],
        ('y', "X") => two("X", c * k, "Z", c * s),
        ('y', "Z") => two("Z", c * k, "X", -c * s),
        ('z', "I") | ('z', "Z") => vec![(p.to_string(), c)],
        ('z', "X") => two("X", c * k, "Y", -c * s),
        ('z', "Y") => two("Y", c * k, "X", c * s),
        _ => unreachable!("unhandled {axis} {p}"),
    }
}

/// Canonicalize a `(word, coeff)` list to a sorted, above-floor set for comparison.
fn canon(mut v: Vec<(String, f64)>) -> Vec<(String, f64)> {
    v.retain(|(_, c)| c.abs() > FLOOR);
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

#[track_caller]
fn assert_terms_close(got: Vec<(String, f64)>, want: Vec<(String, f64)>, ctx: &str) {
    let g = canon(got);
    let w = canon(want);
    assert_eq!(g.len(), w.len(), "{ctx}: term count\ngot={g:?}\nwant={w:?}");
    for (a, b) in g.iter().zip(w.iter()) {
        assert_eq!(a.0, b.0, "{ctx}: word {} vs {}", a.0, b.0);
        assert_close(a.1, b.1, 1e-12);
    }
}

/// A named single-qubit rotation axis paired for the old and new sums.
type AxisGatePair = (char, fn(&mut OldSum, f64), fn(&mut NewSum, f64));

#[test]
fn single_qubit_rotation_exhaustive_vs_reference() {
    // A spread of angles including the commuting/anticommuting extremes.
    let thetas = [0.0, 0.3, 0.7, 1.1, std::f64::consts::FRAC_PI_2, 2.1, -0.9];
    let axes: [AxisGatePair; 3] = [
        ('x', |s, t| s.rx(0, t), |s, t| s.rx(0, t)),
        ('y', |s, t| s.ry(0, t), |s, t| s.ry(0, t)),
        ('z', |s, t| s.rz(0, t), |s, t| s.rz(0, t)),
    ];
    for (axis, gold, gnew) in axes {
        for p in ["I", "X", "Y", "Z"] {
            for &c in &[1.0f64, -1.0, 2.5, -0.75] {
                for &theta in &thetas {
                    let terms = vec![(p.to_string(), c)];
                    let mut old = build_old_sum(1, &terms);
                    let mut new = build_new_sum(1, &terms);
                    gold(&mut old, theta);
                    gnew(&mut new, theta);

                    let want = reference_single(axis, p, c, theta);
                    let ctx = format!("r{axis} {p} c={c} θ={theta}");
                    // NEW matches the independent reference...
                    assert_terms_close(new_support(&new), want.clone(), &ctx);
                    // ...and so does OLD (ties the two crates through the oracle).
                    assert_terms_close(old_support(&old), want, &ctx);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// (c) pauli_error — random probability triples vs old, term for term.
// ---------------------------------------------------------------------------

/// ω-based single-site anticommutation `ω(P, Q)` on the 2-bit codes `(x, z)`.
fn omega(px: bool, pz: bool, qx: bool, qz: bool) -> bool {
    (px & qz) ^ (pz & qx)
}

/// The independent channel eigenvalue `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}` at a site,
/// with `p_I = 1 − pX − pY − pZ`. Pauli codes `(x, z)`: I=(0,0), X=(1,0),
/// Z=(0,1), Y=(1,1).
fn lambda(px: bool, pz: bool, p: [f64; 3]) -> f64 {
    let [pxp, pyp, pzp] = p; // probabilities of X, Y, Z errors
    let pip = 1.0 - pxp - pyp - pzp;
    let paulis = [
        ((false, false), pip), // I
        ((true, false), pxp),  // X
        ((true, true), pyp),   // Y
        ((false, true), pzp),  // Z
    ];
    paulis
        .into_iter()
        .map(
            |((qx, qz), pq)| {
                if omega(px, pz, qx, qz) { -pq } else { pq }
            },
        )
        .sum()
}

#[test]
fn pauli_error_matches_old_and_reference() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            for _ in 0..8 {
                let q = rng.random_range(0..n);
                // A physical probability triple: each in [0, 0.3), so the implicit
                // p_I = 1 − Σ stays positive.
                let p = [
                    rng.random_range(0.0..0.3),
                    rng.random_range(0.0..0.3),
                    rng.random_range(0.0..0.3),
                ];

                let mut old = build_old_sum(n, &terms);
                let mut new = build_new_sum(n, &terms);
                old.pauli_error(q, p);
                new.pauli_error(q, p);

                // Diagonal channel: no re-key, no fan-out — every term keeps its
                // key and is scaled by λ_P. Supports match term for term.
                assert_rot_supports_match(&old, &new, TOL);

                // And the per-term scaling equals the independent ω-based λ_P
                // applied to the original coefficient.
                let before: std::collections::BTreeMap<String, f64> =
                    new_support(&build_new_sum(n, &terms)).into_iter().collect();
                for (word, coeff) in new_support(&new) {
                    let bytes = word.as_bytes();
                    let ch = bytes[q] as char;
                    let (px, pz) = match ch {
                        'I' => (false, false),
                        'X' => (true, false),
                        'Y' => (true, true),
                        'Z' => (false, true),
                        other => panic!("unexpected pauli char {other}"),
                    };
                    let expect = before[&word] * lambda(px, pz, p);
                    assert_close(coeff, expect, TOL.max(expect.abs() * 1e-9));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// (d) COLLISION MERGE — two rotation branches that land on the same key add.
// ---------------------------------------------------------------------------

/// One qubit, two anticommuting terms whose rotation branches cross onto each
/// other, so `accumulate_batch` must *merge* a produced branch onto an existing
/// (cos-scaled) term. Checked against the closed form and against the old crate.
#[test]
fn rotation_branch_collision_merges_correctly() {
    // (axis, (word_a, coeff_a), (word_b, coeff_b),
    //  closed form for the two resulting coefficients as fns of (a, b, sin, cos))
    // rz on {X:b, Y:a}: X ← b·cos + a·sin, Y ← a·cos − b·sin.
    // rx on {Y:a, Z:b}: Y ← a·cos + b·sin, Z ← b·cos − a·sin.
    // ry on {X:b, Z:d}: X ← b·cos − d·sin, Z ← d·cos + b·sin.
    #[allow(clippy::type_complexity)]
    let cases: &[(
        char,
        (&str, f64),
        (&str, f64),
        fn(f64, f64, f64, f64) -> [(&'static str, f64); 2],
    )] = &[
        ('z', ("X", 0.7), ("Y", -0.4), |a, b, s, k| {
            // a=coeff of X (0.7), b=coeff of Y (−0.4)
            [("X", a * k + b * s), ("Y", b * k - a * s)]
        }),
        ('x', ("Y", 0.9), ("Z", -1.3), |a, b, s, k| {
            // a=coeff of Y, b=coeff of Z
            [("Y", a * k + b * s), ("Z", b * k - a * s)]
        }),
        ('y', ("X", -0.6), ("Z", 1.2), |a, b, s, k| {
            // a=coeff of X, b=coeff of Z
            [("X", a * k - b * s), ("Z", b * k + a * s)]
        }),
    ];

    for &theta in &[0.4f64, 1.1, 2.3, -0.8] {
        let (s, k) = (theta.sin(), theta.cos());
        for (axis, (wa, ca), (wb, cb), closed) in cases {
            let terms = vec![(wa.to_string(), *ca), (wb.to_string(), *cb)];
            let mut old = build_old_sum(1, &terms);
            let mut new = build_new_sum(1, &terms);
            match axis {
                'x' => {
                    old.rx(0, theta);
                    new.rx(0, theta);
                }
                'y' => {
                    old.ry(0, theta);
                    new.ry(0, theta);
                }
                _ => {
                    old.rz(0, theta);
                    new.rz(0, theta);
                }
            }

            let want: Vec<(String, f64)> = closed(*ca, *cb, s, k)
                .into_iter()
                .map(|(w, c)| (w.to_string(), c))
                .collect();
            let ctx = format!("collision r{axis} θ={theta}");
            // The merged coefficients match the closed form in the NEW crate...
            assert_terms_close(new_support(&new), want.clone(), &ctx);
            // ...and in the OLD crate.
            assert_terms_close(old_support(&old), want, &ctx);
        }
    }
}

// ---------------------------------------------------------------------------
// (e) THE RESTORED SURFACE — the axis-generic `rotate_1` entry point and the
// `rx_many`/`ry_many`/`rz_many` batch forms, both of which the old
// `ppvm_traits::RotationOne` provided (and the Python bindings call). The new
// trait defaults `rx`/`ry`/`rz` onto `rotate_1` and the batch forms onto those,
// so this pins that the dispatch reaches the same kernel the per-axis fast paths
// use — including the `Pauli::I` axis, which commutes with everything and which
// the old `levi_civita(p, I) = (0, _)` made a no-op.
// ---------------------------------------------------------------------------

#[test]
fn rotate_1_matches_old_on_every_axis() {
    use ppvm_traits::char::Pauli as OldPauli;
    use ppvm_traits_2::Pauli as NewPauli;

    let axes = [
        (OldPauli::I, NewPauli::I),
        (OldPauli::X, NewPauli::X),
        (OldPauli::Y, NewPauli::Y),
        (OldPauli::Z, NewPauli::Z),
    ];
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 8);
            for (old_axis, new_axis) in axes {
                let mut old = build_old_sum(n, &terms);
                let mut new = build_new_sum(n, &terms);
                for q in 0..n {
                    let theta = rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
                    old.rotate_1(old_axis, q, theta);
                    new.rotate_1(new_axis, q, theta);
                    assert_rot_supports_match(&old, &new, TOL);
                }
            }
        }
    }
}

#[test]
fn rotate_1_identity_axis_is_a_no_op() {
    use ppvm_traits_2::Pauli as NewPauli;

    let mut rng = seeded_rng(4242);
    let terms = random_terms(&mut rng, 4, 6);
    let mut new = build_new_sum(4, &terms);
    let before = new_support(&new);
    for q in 0..4 {
        new.rotate_1(NewPauli::I, q, 0.7);
    }
    assert_eq!(new_support(&new), before);
}

#[test]
fn rotation_batch_forms_match_old() {
    // Three rotations *per qubit* fan out 1→2 terms each, so the support grows
    // toward the full `4ⁿ` basis; the widths here keep the per-gate compare
    // bounded while still driving every branch and merge path.
    for &seed in &SEEDS[..3] {
        let mut rng = seeded_rng(seed);
        for &n in WIDTHS.iter().filter(|&&n| n <= 5) {
            let terms = random_terms(&mut rng, n, 8);
            let targets: Vec<usize> = (0..n).rev().collect();
            let theta = rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);

            let mut old = build_old_sum(n, &terms);
            let mut new = build_new_sum(n, &terms);
            old.rx_many(&targets, theta);
            new.rx_many(&targets, theta);
            assert_rot_supports_match(&old, &new, TOL);

            old.ry_many(&targets, theta);
            new.ry_many(&targets, theta);
            assert_rot_supports_match(&old, &new, TOL);

            old.rz_many(&targets, theta);
            new.rz_many(&targets, theta);
            assert_rot_supports_match(&old, &new, TOL);

            // The batch form is exactly the per-target loop, in order.
            let mut loops = build_new_sum(n, &terms);
            for &q in &targets {
                loops.rx(q, theta);
            }
            for &q in &targets {
                loops.ry(q, theta);
            }
            for &q in &targets {
                loops.rz(q, theta);
            }
            assert_eq!(new_support(&loops), new_support(&new));
        }
    }
}
