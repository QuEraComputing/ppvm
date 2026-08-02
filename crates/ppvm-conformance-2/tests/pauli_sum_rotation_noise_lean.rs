// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for the **non-Clifford** surface of
//! `ppvm-pauli-sum-2::PauliSum<f64>`: the rotation branch of
//! `lean/PPVM/Instantiations/Rotation.lean` and the unital Pauli-channel
//! eigenvalue of `lean/PPVM/Algebra/Noise.lean`, reproduced as randomized (and,
//! where the single-qubit case is finite, exhaustive) property tests on the sum.
//!
//! Laws reproduced:
//! * `anticommute_new_key` — when `{G, P} = 0` the branch key `iGP` is a
//!   *genuinely new* Pauli with symplectic bits `G ⊕ P`, distinct from both `G`
//!   and `P`; exhaustive over the single-qubit `(G, P)` grid.
//! * `commute_bits` / `rot_zero` — the commuting case is inert (no new key), and
//!   `θ = 0` is the identity rotation.
//! * `rot_norm_sq` — a rotation preserves the ℓ² norm of the coefficient vector
//!   (it is an orthogonal map on the anticommuting `(P, iGP)` plane).
//! * `rot_neg_rot` — reversibility: `R_{−θ} ∘ R_θ = id`.
//! * `rot_rot` — angle addition: `R_θ ∘ R_φ = R_{θ+φ}`; in particular a Trotter
//!   step of two `rz(θ)` equals one `rz(2θ)`.
//! * `pauli_channel_eigenvalue_omega` — the diagonal channel scales `P` by
//!   `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)} = 1 − 2·Σ_{Q anti} p_Q`, tied to the symplectic
//!   form `ω`; exhaustive over the single-qubit `P`, randomized over `p`.

use ppvm_conformance_2::{
    NewKey, NewSum, assert_close, build_new_sum, new_support, random_terms, seeded_rng,
};
use ppvm_traits_2::{PauliError, RotationOne, Word};

use rand::RngExt;
use rand::rngs::StdRng;
use std::collections::BTreeMap;

const SEEDS: [u64; 10] = [1, 2, 3, 7, 42, 99, 123, 777, 2024, 31337];
const WIDTHS: [usize; 6] = [1, 2, 3, 5, 8, 12];
const TOL: f64 = 1e-9;
/// Coefficients below this magnitude are treated as *absent*. A rotation branch
/// that cancels back leaves a physical near-zero residue on the branch key (the
/// `reduce` step drops only *exact* zeros), so composition/inverse identities are
/// compared on the above-floor support. Well below any `O(1)` coefficient.
const FLOOR: f64 = 1e-9;

/// The sum's support as a `key → coeff` map (canonical strings).
fn support_map(s: &NewSum) -> BTreeMap<String, f64> {
    new_support(s).into_iter().collect()
}

/// The above-floor support as a `key → coeff` map (drops branch residue).
fn above_floor_map(s: &NewSum) -> BTreeMap<String, f64> {
    new_support(s)
        .into_iter()
        .filter(|(_, c)| c.abs() > FLOOR)
        .collect()
}

/// The squared ℓ² norm of the coefficient vector.
fn norm_sq(s: &NewSum) -> f64 {
    new_support(s).into_iter().map(|(_, c)| c * c).sum()
}

/// A random rotation angle in `(−π, π)`.
fn angle(rng: &mut StdRng) -> f64 {
    rng.random_range(-std::f64::consts::PI..std::f64::consts::PI)
}

/// Apply `r{axis}` (`axis` ∈ `{0:x, 1:y, 2:z}`) on `q` by `theta`.
fn apply_axis(s: &mut NewSum, axis: usize, q: usize, theta: f64) {
    match axis {
        0 => s.rx(q, theta),
        1 => s.ry(q, theta),
        _ => s.rz(q, theta),
    }
}

// ---------------------------------------------------------------------------
// rot_norm_sq — a rotation preserves the ℓ² norm of the coefficient vector.
// ---------------------------------------------------------------------------

#[test]
fn rotation_preserves_l2_norm() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 25);
            let mut s = build_new_sum(n, &terms);
            let before = norm_sq(&s);
            // A run of rotations about random axes/qubits/angles; each is an
            // orthogonal map on the coefficient space, so ‖c‖² is invariant. The
            // support fans out toward `4ⁿ`, so cap the run on wide registers.
            let rounds = if n <= 5 { 25 } else { 12 };
            for _ in 0..rounds {
                let axis = rng.random_range(0..3usize);
                let q = rng.random_range(0..n);
                apply_axis(&mut s, axis, q, angle(&mut rng));
                assert_close(norm_sq(&s), before, 1e-7 * (1.0 + before));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// rot_neg_rot — reversibility: R_{−θ} ∘ R_θ = id.
// ---------------------------------------------------------------------------

#[test]
fn rotation_is_reversible() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);
            let base = build_new_sum(n, &terms);
            let base_map = above_floor_map(&base);

            for axis in 0..3usize {
                let q = rng.random_range(0..n);
                let theta = angle(&mut rng);
                let mut s = build_new_sum(n, &terms);
                apply_axis(&mut s, axis, q, theta);
                apply_axis(&mut s, axis, q, -theta);
                // Back to the original support, key by key.
                let m = above_floor_map(&s);
                assert_eq!(
                    m.len(),
                    base_map.len(),
                    "R_-θ∘R_θ changed support size (axis {axis} q {q} θ {theta})"
                );
                for (k, v) in &base_map {
                    let got = m.get(k).unwrap_or_else(|| panic!("key {k} vanished"));
                    assert_close(*got, *v, 1e-8);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// rot_rot — angle addition: R_θ ∘ R_φ = R_{θ+φ}; Trotter two rz(θ) = one rz(2θ).
// ---------------------------------------------------------------------------

#[test]
fn rotation_angle_addition_and_trotter() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &WIDTHS {
            let terms = random_terms(&mut rng, n, 30);

            for axis in 0..3usize {
                let q = rng.random_range(0..n);
                let theta = angle(&mut rng);
                let phi = angle(&mut rng);

                // R_θ ∘ R_φ.
                let mut composed = build_new_sum(n, &terms);
                apply_axis(&mut composed, axis, q, phi);
                apply_axis(&mut composed, axis, q, theta);

                // R_{θ+φ}. Keep θ+φ in a range where the branch stays regular.
                let mut single = build_new_sum(n, &terms);
                apply_axis(&mut single, axis, q, theta + phi);

                let a = above_floor_map(&composed);
                let b = above_floor_map(&single);
                assert_eq!(
                    a.len(),
                    b.len(),
                    "angle-addition support size mismatch (axis {axis} q {q})"
                );
                for (k, va) in &a {
                    let vb = b.get(k).unwrap_or_else(|| panic!("key {k} missing"));
                    assert_close(*va, *vb, 1e-8);
                }
            }

            // The Trotter special case, spelled out: two rz(θ) == one rz(2θ).
            let q = rng.random_range(0..n);
            let theta = angle(&mut rng) * 0.5; // so 2θ stays in (−π, π)
            let mut two = build_new_sum(n, &terms);
            two.rz(q, theta);
            two.rz(q, theta);
            let mut one = build_new_sum(n, &terms);
            one.rz(q, 2.0 * theta);
            let a = above_floor_map(&two);
            let b = above_floor_map(&one);
            assert_eq!(a.len(), b.len(), "Trotter rz support size mismatch");
            for (k, va) in &a {
                assert_close(*va, b[k], 1e-8);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// anticommute_new_key / commute_bits / rot_zero — exhaustive single-qubit.
// ---------------------------------------------------------------------------

/// Symplectic bits `(x, z)` of a single-qubit Pauli char.
fn bits_of(ch: char) -> (bool, bool) {
    match ch {
        'I' => (false, false),
        'X' => (true, false),
        'Y' => (true, true),
        'Z' => (false, true),
        _ => panic!("not a pauli: {ch}"),
    }
}

/// The rotation-axis generator `G` bits per axis (`0:x, 1:y, 2:z`).
fn axis_bits(axis: usize) -> (bool, bool) {
    match axis {
        0 => (true, false), // X
        1 => (true, true),  // Y
        _ => (false, true), // Z
    }
}

#[test]
fn anticommuting_branch_produces_distinct_new_key_with_xor_bits() {
    // θ = π/2: both cos and sin are nonzero, so an anticommuting term yields two
    // terms (the original `P` and the branch `iGP`); a commuting term stays one.
    let theta = std::f64::consts::FRAC_PI_2;
    for axis in 0..3usize {
        let (gx, gz) = axis_bits(axis);
        for p in ["I", "X", "Y", "Z"] {
            let (px, pz) = bits_of(p.chars().next().unwrap());
            let anti = (gx & pz) ^ (gz & px); // ω(G, P)

            let terms = vec![(p.to_string(), 1.0)];
            let mut s = build_new_sum(1, &terms);
            apply_axis(&mut s, axis, 0, theta);
            let support = new_support(&s);

            if anti {
                // Exactly two terms: P (cos) and the branch key (sin).
                assert_eq!(
                    support.len(),
                    2,
                    "axis {axis} P {p}: anticommuting branch should fan out to 2"
                );
                let branch: Vec<&(String, f64)> = support.iter().filter(|(w, _)| w != p).collect();
                assert_eq!(branch.len(), 1, "axis {axis} P {p}: one new key");
                let new_key = branch[0].0.chars().next().unwrap();
                let (nx, nz) = bits_of(new_key);
                // Bits are exactly G ⊕ P (the mulBits key product).
                assert_eq!((nx, nz), (gx ^ px, gz ^ pz), "axis {axis} P {p}: XOR bits");
                // ...and the new key differs from both G and P.
                let g_char = match (gx, gz) {
                    (true, false) => 'X',
                    (true, true) => 'Y',
                    _ => 'Z',
                };
                assert_ne!(new_key, p.chars().next().unwrap(), "new key == P");
                assert_ne!(new_key, g_char, "new key == G");
            } else {
                // Commuting: rot is inert on the key set — the single P survives.
                assert_eq!(
                    support.len(),
                    1,
                    "axis {axis} P {p}: commuting term must not branch"
                );
                assert_eq!(support[0].0, p, "axis {axis} P {p}: key unchanged");
            }
        }
    }
}

#[test]
fn rotation_by_zero_is_identity() {
    // rot_zero: R_0 = id on any support and any axis.
    //
    // "Identity" is identity *as an element of `C[K]`*, i.e. on the reduced
    // support. At θ = 0 the branch coefficient is `sinθ · c = 0`, and the engine
    // merges that zero branch rather than skipping it — deliberately, because old
    // does (`map_insert` pass 2 calls `add_assign` for every produced term, and
    // `add_assign` inserts a 0.0). So `R_0` can leave an extra *zero-coefficient*
    // key, which `reduce` removes; both sides are compared in reduced form, and
    // `reduce_structural` (`lean/PPVM/Algebra/GradedMap.lean`) is what says the two
    // representations denote the same map.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 3, 8] {
            let terms = random_terms(&mut rng, n, 20);
            let mut base_sum = build_new_sum(n, &terms);
            base_sum.reduce();
            let base = support_map(&base_sum);
            for axis in 0..3usize {
                let q = rng.random_range(0..n);
                let mut s = build_new_sum(n, &terms);
                apply_axis(&mut s, axis, q, 0.0);
                s.reduce();
                assert_eq!(support_map(&s), base, "R_0 not identity (axis {axis})");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// pauli_channel_eigenvalue_omega — λ_P = Σ_Q p_Q (−1)^{ω(P,Q)} = 1 − 2Σ_anti p_Q.
// ---------------------------------------------------------------------------

/// The raw eigenvalue `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}` with `p_I = 1 − Σ`.
fn lambda_raw(px: bool, pz: bool, p: [f64; 3]) -> f64 {
    let [pxp, pyp, pzp] = p;
    let pip = 1.0 - pxp - pyp - pzp;
    let qs = [
        ((false, false), pip),
        ((true, false), pxp),
        ((true, true), pyp),
        ((false, true), pzp),
    ];
    qs.into_iter()
        .map(|((qx, qz), pq)| {
            let anti = (px & qz) ^ (pz & qx);
            if anti { -pq } else { pq }
        })
        .sum()
}

/// The collapsed eigenvalue `λ_P = 1 − 2·Σ_{Q anti} p_Q`.
fn lambda_collapsed(px: bool, pz: bool, p: [f64; 3]) -> f64 {
    let [pxp, pyp, pzp] = p;
    let pip = 1.0 - pxp - pyp - pzp;
    let qs = [
        ((false, false), pip),
        ((true, false), pxp),
        ((true, true), pyp),
        ((false, true), pzp),
    ];
    let anti_sum: f64 = qs
        .into_iter()
        .filter(|&((qx, qz), _)| (px & qz) ^ (pz & qx))
        .map(|(_, pq)| pq)
        .sum();
    1.0 - 2.0 * anti_sum
}

#[test]
fn pauli_channel_eigenvalue_omega_matches() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for _ in 0..64 {
            let p = [
                rng.random_range(0.0..0.3),
                rng.random_range(0.0..0.3),
                rng.random_range(0.0..0.3),
            ];
            // Exhaustive over the single-qubit Pauli P.
            for pstr in ["I", "X", "Y", "Z"] {
                let (px, pz) = bits_of(pstr.chars().next().unwrap());

                // The two Lean forms of λ_P agree (this IS the theorem's identity).
                let raw = lambda_raw(px, pz, p);
                let collapsed = lambda_collapsed(px, pz, p);
                assert_close(raw, collapsed, 1e-12);

                // The channel scales a pure P by exactly λ_P.
                let mut s = build_new_sum(1, &[(pstr.to_string(), 1.0)]);
                s.pauli_error(0, p);
                let got = s
                    .get(&NewKey::from(pstr))
                    .expect("diagonal channel keeps the key");
                assert_close(got, collapsed, TOL);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// A `Word`-trait smoke over the branch: the produced keys are valid n-site words.
// ---------------------------------------------------------------------------

#[test]
fn branch_keys_have_correct_arity() {
    let mut rng = seeded_rng(7);
    for &n in &[2usize, 4, 8] {
        let terms = random_terms(&mut rng, n, 10);
        let mut s = build_new_sum(n, &terms);
        s.rx(0, 1.2);
        s.ry(1 % n, 0.6);
        for (w, _) in new_support(&s) {
            let key = NewKey::from(w.as_str());
            assert_eq!(key.n_sites(), n, "branch key {w} lost its arity");
        }
    }
}
