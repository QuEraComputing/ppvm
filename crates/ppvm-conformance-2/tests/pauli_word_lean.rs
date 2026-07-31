// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for `ppvm-pauli-word-2::PauliWord`: the
//! machine-checked semantics of `lean/PPVM/Pauli/**` reproduced as Rust tests,
//! grounded in a genuine ℤ[i] (Gaussian-integer) matrix reference exactly as the
//! Lean development grounds `phaseExp` in `PPVM.PauliMatrix.pauliMat_mul`.
//!
//! Coverage (exhaustive for finite single-qubit cases, randomized for n-qubit):
//!
//! * `Phase.lean` `phaseExp_eq_ref` / `Matrix.lean` `pauliMat_mul` — the new
//!   `key_mul` phase exponent equals the base-`i` exponent of the real 2×2 (and
//!   n-fold Kronecker) ℤ[i] matrix product.
//! * `Word.lean` `phaseExpN_cocycle` — the summed phase is a 2-cocycle
//!   (associativity of the twisted product).
//! * `Word.lean` `phaseExpN_self` — `P·P = +I`.
//! * `Word.lean` `phaseExpN_sub_comm` — `P·Q = (−1)^{ω} Q·P`.
//! * `Conjugation.lean` `conjH_*`/`conjS_*`/`conjCNOT_*`/`conjCZ_*` — the
//!   Clifford conjugation tables (HXH=Z, HYH=−Y, SXS†=Y, and the CNOT/CZ
//!   generator tables), reproduced as the Lean `conj*` functions, checked
//!   against those tables, shown to be group homomorphisms, and grounded in
//!   genuine ℤ[i] gate-conjugation matrices; the *bit* half is tied to the crate
//!   under test's Clifford map.
//! * `Symplectic.lean` `hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/
//!   `czAct_isometry` — every generator preserves the symplectic form `ω`, run
//!   on the crate's own `Clifford` implementation.

use ppvm_conformance_2::{random_pauli_string, seeded_rng};
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{Clifford, KeyProduct, PauliBits, Phase, Word};
use rand::RngExt;

type New = NewWord<u64>;

const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];

// ===========================================================================
// ℤ[i] matrix reference (mirrors lean/PPVM/Pauli/Matrix.lean)
// ===========================================================================

/// A Gaussian integer `ℤ[i]` (exact; no floats).
type Z = num::Complex<i64>;
/// A dense matrix over `ℤ[i]`.
type Mat = Vec<Vec<Z>>;

#[inline]
fn z(re: i64, im: i64) -> Z {
    num::Complex::new(re, im)
}

/// `iᵏ ∈ ℤ[i]` for exponent `k` — the scalar the twisted product emits.
fn ipow(k: u8) -> Z {
    match k & 3 {
        0 => z(1, 0),
        1 => z(0, 1),
        2 => z(-1, 0),
        _ => z(0, -1),
    }
}

/// The single-qubit Pauli `g(x,z) = iˣᶻ Xˣ Zᶻ` as a genuine 2×2 ℤ[i] matrix,
/// with `Y = iXZ` (matches `PauliMatrix.pauliMat`).
fn pauli_mat(x: bool, zbit: bool) -> Mat {
    match (x, zbit) {
        (false, false) => vec![vec![z(1, 0), z(0, 0)], vec![z(0, 0), z(1, 0)]],
        (true, false) => vec![vec![z(0, 0), z(1, 0)], vec![z(1, 0), z(0, 0)]],
        (false, true) => vec![vec![z(1, 0), z(0, 0)], vec![z(0, 0), z(-1, 0)]],
        (true, true) => vec![vec![z(0, 0), z(0, -1)], vec![z(0, 1), z(0, 0)]],
    }
}

fn pauli_mat_of_char(c: char) -> Mat {
    match c {
        'I' => pauli_mat(false, false),
        'X' => pauli_mat(true, false),
        'Z' => pauli_mat(false, true),
        'Y' => pauli_mat(true, true),
        other => panic!("bad Pauli char {other}"),
    }
}

fn matmul(a: &Mat, b: &Mat) -> Mat {
    let (n, m, p) = (a.len(), b.len(), b[0].len());
    let mut out = vec![vec![z(0, 0); p]; n];
    for (i, orow) in out.iter_mut().enumerate() {
        for (k, brow) in b.iter().enumerate().take(m) {
            let aik = a[i][k];
            for (j, oj) in orow.iter_mut().enumerate().take(p) {
                *oj += aik * brow[j];
            }
        }
    }
    out
}

/// Kronecker product `a ⊗ b` (qubit-`0`-outermost convention).
fn kron(a: &Mat, b: &Mat) -> Mat {
    let (ar, ac, br, bc) = (a.len(), a[0].len(), b.len(), b[0].len());
    let mut out = vec![vec![z(0, 0); ac * bc]; ar * br];
    for i in 0..ar {
        for j in 0..ac {
            for k in 0..br {
                for l in 0..bc {
                    out[i * br + k][j * bc + l] = a[i][j] * b[k][l];
                }
            }
        }
    }
    out
}

/// Full n-fold Kronecker matrix of a Pauli string (qubit `0` leftmost).
fn word_mat(s: &str) -> Mat {
    let mut acc: Mat = vec![vec![z(1, 0)]];
    for c in s.chars() {
        acc = kron(&acc, &pauli_mat_of_char(c));
    }
    acc
}

fn scalar_mul(s: Z, a: &Mat) -> Mat {
    a.iter()
        .map(|row| row.iter().map(|&e| s * e).collect())
        .collect()
}

// ===========================================================================
// phaseExp_eq_ref / pauliMat_mul
// ===========================================================================

const LETTERS: [char; 4] = ['I', 'X', 'Y', 'Z'];

#[test]
fn phase_exp_equals_matrix_exponent_single_qubit_exhaustive() {
    // Matrix.lean pauliMat_mul: g(p)·g(q) = i^{phaseExp} · g(p⊕q), all 16 cases.
    for p in LETTERS {
        for q in LETTERS {
            let a: New = p.to_string().as_str().into();
            let b: New = q.to_string().as_str().into();
            let (r, phase): (New, Phase) = a.key_mul(&b);

            let lhs = matmul(&pauli_mat_of_char(p), &pauli_mat_of_char(q));
            let rhs = scalar_mul(ipow(phase.exponent()), &word_mat(&r.to_string()));
            assert_eq!(
                lhs, rhs,
                "{p}·{q}: key_mul phase disagrees with ℤ[i] matrix"
            );
        }
    }
}

#[test]
fn phase_exp_equals_matrix_exponent_n_qubit_random() {
    // Same identity for n-fold Kronecker matrices (n ≤ 5 keeps 2ⁿ ≤ 32).
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=5usize {
            let v = random_pauli_string(&mut rng, n);
            let w = random_pauli_string(&mut rng, n);
            let a: New = v.as_str().into();
            let b: New = w.as_str().into();
            let (r, phase): (New, Phase) = a.key_mul(&b);

            let lhs = matmul(&word_mat(&v), &word_mat(&w));
            let rhs = scalar_mul(ipow(phase.exponent()), &word_mat(&r.to_string()));
            assert_eq!(lhs, rhs, "seed {seed}: {v}·{w} disagrees with ℤ[i] matrix");
        }
    }
}

// ===========================================================================
// phaseExpN_cocycle / phaseExpN_self / phaseExpN_sub_comm
// ===========================================================================

#[test]
fn phase_cocycle_associativity_random() {
    // Word.lean phaseExpN_cocycle: (u·v)·w and u·(v·w) agree on bits AND folded
    // phase — the 2-cocycle identity that makes the twisted product associative.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let u: New = random_pauli_string(&mut rng, n).as_str().into();
            let v: New = random_pauli_string(&mut rng, n).as_str().into();
            let w: New = random_pauli_string(&mut rng, n).as_str().into();

            let (uv, p1) = u.key_mul(&v);
            let (uv_w, p2) = uv.key_mul(&w);
            let left_phase = p1.compose(p2);

            let (vw, p3) = v.key_mul(&w);
            let (u_vw, p4) = u.key_mul(&vw);
            let right_phase = p3.compose(p4);

            assert_eq!(uv_w, u_vw, "seed {seed} n {n}: associativity bits");
            assert_eq!(left_phase, right_phase, "seed {seed} n {n}: cocycle phase");
        }
    }
}

#[test]
fn square_is_plus_identity_random() {
    // Word.lean phaseExpN_self: P·P = +I.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let s = random_pauli_string(&mut rng, n);
            let p: New = s.as_str().into();
            let (r, phase) = p.key_mul(&p);
            assert_eq!(r, New::new(n), "seed {seed}: {s}² word != I");
            assert_eq!(phase, Phase::Pos1, "seed {seed}: {s}² phase != +1");
        }
    }
}

/// Symplectic form `ω(p,q) = Σᵢ (x_pᵢ·z_qᵢ ⊕ z_pᵢ·x_qᵢ)` over 𝔽₂, read from the
/// crate's own bit accessors (matches `Symplectic.lean omegaFun`).
fn omega(p: &New, q: &New) -> u8 {
    let mut acc = 0u8;
    for i in 0..p.n_sites() {
        acc ^= (p.x_bit(i) && q.z_bit(i)) as u8;
        acc ^= (p.z_bit(i) && q.x_bit(i)) as u8;
    }
    acc
}

#[test]
fn commutation_is_symplectic_form_random() {
    // Word.lean phaseExpN_sub_comm: P·Q·(Q·P)⁻¹ = (−1)^{ω(P,Q)}.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let p: New = random_pauli_string(&mut rng, n).as_str().into();
            let q: New = random_pauli_string(&mut rng, n).as_str().into();
            let (rpq, pq) = p.key_mul(&q);
            let (rqp, qp) = q.key_mul(&p);
            assert_eq!(rpq, rqp, "seed {seed}: P·Q and Q·P differ on bits");

            let expected = if omega(&p, &q) == 0 {
                Phase::Pos1
            } else {
                Phase::Neg1
            };
            assert_eq!(
                pq.compose(qp.inverse()),
                expected,
                "seed {seed}: comm phase"
            );
        }
    }
}

// ===========================================================================
// Clifford conjugation tables (mirrors lean/PPVM/Pauli/Conjugation.lean)
// ===========================================================================

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct P1 {
    phase: u8,
    x: bool,
    z: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct P2 {
    phase: u8,
    xc: bool,
    zc: bool,
    xt: bool,
    zt: bool,
}

/// Single-qubit product phase exponent (the `phase/mul.rs` booleans, = Lean
/// `phaseExp`). Kept as the verbatim `sign`/`imag` sum-of-products the kernel and
/// `PPVM.PauliPhase.signBit`/`imagBit` use, so `clippy::nonminimal_bool` is
/// silenced to preserve the one-to-one correspondence with the spec.
#[allow(clippy::nonminimal_bool)]
fn phase_exp(a: bool, b: bool, c: bool, d: bool) -> u8 {
    let sign = (a && b && c && !d) || (a && !b && !c && d) || (!a && b && c && d);
    let imag = (a && !b && d) || (a && !c && d) || (!a && b && c) || (b && c && !d);
    (2 * sign as u8 + imag as u8) % 4
}

fn mul_p1(p: P1, q: P1) -> P1 {
    P1 {
        phase: (p.phase + q.phase + phase_exp(p.x, p.z, q.x, q.z)) % 4,
        x: p.x ^ q.x,
        z: p.z ^ q.z,
    }
}

fn mul_p2(p: P2, q: P2) -> P2 {
    P2 {
        phase: (p.phase
            + q.phase
            + phase_exp(p.xc, p.zc, q.xc, q.zc)
            + phase_exp(p.xt, p.zt, q.xt, q.zt))
            % 4,
        xc: p.xc ^ q.xc,
        zc: p.zc ^ q.zc,
        xt: p.xt ^ q.xt,
        zt: p.zt ^ q.zt,
    }
}

/// `conjH` (Conjugation.lean): swap `(x,z)`, sign `2·(x∧z)`.
fn conj_h(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x && p.z { 2 } else { 0 }) % 4,
        x: p.z,
        z: p.x,
    }
}

/// `conjS` (Conjugation.lean): `(x,z) ↦ (x, x⊕z)`, sign `2·(x∧z)`.
fn conj_s(p: P1) -> P1 {
    P1 {
        phase: (p.phase + if p.x && p.z { 2 } else { 0 }) % 4,
        x: p.x,
        z: p.x ^ p.z,
    }
}

/// `conjCNOT` (Conjugation.lean): `x_t ⊕= x_c`, `z_c ⊕= z_t`, sign `2` iff
/// `x_c ∧ z_t ∧ (x_t = z_c)`.
fn conj_cnot(p: P2) -> P2 {
    let delta = if p.xc && p.zt && (p.xt == p.zc) { 2 } else { 0 };
    P2 {
        phase: (p.phase + delta) % 4,
        xc: p.xc,
        zc: p.zc ^ p.zt,
        xt: p.xt ^ p.xc,
        zt: p.zt,
    }
}

/// `conjCZ` (Conjugation.lean): `z_c ⊕= x_t`, `z_t ⊕= x_c`, sign `2` iff
/// `x_c ∧ x_t ∧ (z_c ≠ z_t)`.
fn conj_cz(p: P2) -> P2 {
    let delta = if p.xc && p.xt && (p.zc ^ p.zt) { 2 } else { 0 };
    P2 {
        phase: (p.phase + delta) % 4,
        xc: p.xc,
        zc: p.zc ^ p.xt,
        xt: p.xt,
        zt: p.zt ^ p.xc,
    }
}

#[test]
fn conjugation_generator_tables() {
    // Exact `by decide` theorems of Conjugation.lean, as assertions.
    // H: HXH=Z, HZH=X, HYH=−Y.
    assert_eq!(
        conj_h(P1 {
            phase: 0,
            x: true,
            z: false
        }),
        P1 {
            phase: 0,
            x: false,
            z: true
        }
    );
    assert_eq!(
        conj_h(P1 {
            phase: 0,
            x: false,
            z: true
        }),
        P1 {
            phase: 0,
            x: true,
            z: false
        }
    );
    assert_eq!(
        conj_h(P1 {
            phase: 0,
            x: true,
            z: true
        }),
        P1 {
            phase: 2,
            x: true,
            z: true
        }
    );
    // S: SXS†=Y, SZS†=Z, SYS†=−X.
    assert_eq!(
        conj_s(P1 {
            phase: 0,
            x: true,
            z: false
        }),
        P1 {
            phase: 0,
            x: true,
            z: true
        }
    );
    assert_eq!(
        conj_s(P1 {
            phase: 0,
            x: false,
            z: true
        }),
        P1 {
            phase: 0,
            x: false,
            z: true
        }
    );
    assert_eq!(
        conj_s(P1 {
            phase: 0,
            x: true,
            z: true
        }),
        P1 {
            phase: 2,
            x: true,
            z: false
        }
    );

    // CNOT generator table: X_c→X_cX_t, X_t→X_t, Z_c→Z_c, Z_t→Z_cZ_t, Y_cY_t→(−1)…
    let g = |xc, zc, xt, zt| P2 {
        phase: 0,
        xc,
        zc,
        xt,
        zt,
    };
    assert_eq!(
        conj_cnot(g(true, false, false, false)),
        g(true, false, true, false)
    );
    assert_eq!(
        conj_cnot(g(false, false, true, false)),
        g(false, false, true, false)
    );
    assert_eq!(
        conj_cnot(g(false, true, false, false)),
        g(false, true, false, false)
    );
    assert_eq!(
        conj_cnot(g(false, false, false, true)),
        g(false, true, false, true)
    );
    assert_eq!(
        conj_cnot(g(true, true, true, true)),
        P2 {
            phase: 2,
            xc: true,
            zc: false,
            xt: false,
            zt: true
        }
    );
    // CZ generator table: X_c→X_cZ_t, X_t→Z_cX_t, Z_c→Z_c, Z_t→Z_t, Y_cX_t→(−1)…
    assert_eq!(
        conj_cz(g(true, false, false, false)),
        g(true, false, false, true)
    );
    assert_eq!(
        conj_cz(g(false, false, true, false)),
        g(false, true, true, false)
    );
    assert_eq!(
        conj_cz(g(false, true, false, false)),
        g(false, true, false, false)
    );
    assert_eq!(
        conj_cz(g(false, false, false, true)),
        g(false, false, false, true)
    );
    assert_eq!(
        conj_cz(g(true, true, true, false)),
        P2 {
            phase: 2,
            xc: true,
            zc: false,
            xt: true,
            zt: true
        }
    );
}

#[test]
fn conjugation_maps_are_group_homomorphisms() {
    // conjH_mul_raw / conjS_mul_raw (16²) and conjCNOT/conjCZ_mul_raw (64²),
    // exhaustive over the bit patterns (phase = 0 on operands; deltas add).
    let p1s: Vec<P1> = (0..4)
        .map(|b| P1 {
            phase: 0,
            x: b & 1 != 0,
            z: b & 2 != 0,
        })
        .collect();
    for &p in &p1s {
        for &q in &p1s {
            assert_eq!(
                conj_h(mul_p1(p, q)),
                mul_p1(conj_h(p), conj_h(q)),
                "conjH hom"
            );
            assert_eq!(
                conj_s(mul_p1(p, q)),
                mul_p1(conj_s(p), conj_s(q)),
                "conjS hom"
            );
        }
    }

    let p2s: Vec<P2> = (0..16)
        .map(|b| P2 {
            phase: 0,
            xc: b & 1 != 0,
            zc: b & 2 != 0,
            xt: b & 4 != 0,
            zt: b & 8 != 0,
        })
        .collect();
    for &p in &p2s {
        for &q in &p2s {
            assert_eq!(
                conj_cnot(mul_p2(p, q)),
                mul_p2(conj_cnot(p), conj_cnot(q)),
                "cnot hom"
            );
            assert_eq!(
                conj_cz(mul_p2(p, q)),
                mul_p2(conj_cz(p), conj_cz(q)),
                "cz hom"
            );
        }
    }
}

// The exact Clifford gate matrices over ℤ[i] (H is √2-scaled to stay integral).
fn sqrt2_h() -> Mat {
    vec![vec![z(1, 0), z(1, 0)], vec![z(1, 0), z(-1, 0)]]
}
fn s_gate() -> Mat {
    vec![vec![z(1, 0), z(0, 0)], vec![z(0, 0), z(0, 1)]]
}
fn s_dag() -> Mat {
    vec![vec![z(1, 0), z(0, 0)], vec![z(0, 0), z(0, -1)]]
}
fn cnot_gate() -> Mat {
    vec![
        vec![z(1, 0), z(0, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(1, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(0, 0), z(0, 0), z(1, 0)],
        vec![z(0, 0), z(0, 0), z(1, 0), z(0, 0)],
    ]
}
fn cz_gate() -> Mat {
    vec![
        vec![z(1, 0), z(0, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(1, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(0, 0), z(1, 0), z(0, 0)],
        vec![z(0, 0), z(0, 0), z(0, 0), z(-1, 0)],
    ]
}

#[test]
fn conjugation_reference_is_grounded_in_zi_matrices() {
    // Ground the Lean `conj*` functions in genuine G·P·G† matrices over ℤ[i],
    // exactly as Matrix.lean grounds phaseExp — closing the "hand-derived
    // formula" caveat for the conjugation signs too.
    let p1s = [
        P1 {
            phase: 0,
            x: false,
            z: false,
        },
        P1 {
            phase: 0,
            x: true,
            z: false,
        },
        P1 {
            phase: 0,
            x: false,
            z: true,
        },
        P1 {
            phase: 0,
            x: true,
            z: true,
        },
    ];
    for &p in &p1s {
        let mp = pauli_mat(p.x, p.z);
        // H: (√2 H)·P·(√2 H) = 2·H·P·H (H Hermitian); expect 2·i^{phase}·result.
        let r = conj_h(p);
        let lhs = matmul(&matmul(&sqrt2_h(), &mp), &sqrt2_h());
        let rhs = scalar_mul(z(2, 0) * ipow(r.phase), &pauli_mat(r.x, r.z));
        assert_eq!(lhs, rhs, "conjH ℤ[i] grounding for {p:?}");
        // S: S·P·S†.
        let r = conj_s(p);
        let lhs = matmul(&matmul(&s_gate(), &mp), &s_dag());
        let rhs = scalar_mul(ipow(r.phase), &pauli_mat(r.x, r.z));
        assert_eq!(lhs, rhs, "conjS ℤ[i] grounding for {p:?}");
    }

    // Two-qubit: CNOT/CZ are Hermitian, self-inverse, integral → no scaling.
    for b in 0..16u8 {
        let p = P2 {
            phase: 0,
            xc: b & 1 != 0,
            zc: b & 2 != 0,
            xt: b & 4 != 0,
            zt: b & 8 != 0,
        };
        let mp = kron(&pauli_mat(p.xc, p.zc), &pauli_mat(p.xt, p.zt));

        let r = conj_cnot(p);
        let lhs = matmul(&matmul(&cnot_gate(), &mp), &cnot_gate());
        let rr = kron(&pauli_mat(r.xc, r.zc), &pauli_mat(r.xt, r.zt));
        assert_eq!(
            lhs,
            scalar_mul(ipow(r.phase), &rr),
            "conjCNOT ℤ[i] grounding {p:?}"
        );

        let r = conj_cz(p);
        let lhs = matmul(&matmul(&cz_gate(), &mp), &cz_gate());
        let rr = kron(&pauli_mat(r.xc, r.zc), &pauli_mat(r.xt, r.zt));
        assert_eq!(
            lhs,
            scalar_mul(ipow(r.phase), &rr),
            "conjCZ ℤ[i] grounding {p:?}"
        );
    }
}

#[test]
fn crate_clifford_bit_map_matches_conjugation_oracle() {
    // Tie the crate under test to the oracle's BIT half: the bare word's
    // Clifford realizes exactly `conj*`'s symplectic action (the phase half is a
    // documented no-op on a phaseless word).
    for b in 0..4u8 {
        let (x, zb) = (b & 1 != 0, b & 2 != 0);
        let s = match (x, zb) {
            (false, false) => "I",
            (true, false) => "X",
            (false, true) => "Z",
            (true, true) => "Y",
        };
        // H.
        let mut w: New = s.into();
        w.h(0);
        let r = conj_h(P1 { phase: 0, x, z: zb });
        assert_eq!((w.x_bit(0), w.z_bit(0)), (r.x, r.z), "H bit map {s}");
        // S.
        let mut w: New = s.into();
        w.s(0);
        let r = conj_s(P1 { phase: 0, x, z: zb });
        assert_eq!((w.x_bit(0), w.z_bit(0)), (r.x, r.z), "S bit map {s}");
    }

    for b in 0..16u8 {
        let (xc, zc, xt, zt) = (b & 1 != 0, b & 2 != 0, b & 4 != 0, b & 8 != 0);
        let letter = |x, zbit| match (x, zbit) {
            (false, false) => 'I',
            (true, false) => 'X',
            (false, true) => 'Z',
            (true, true) => 'Y',
        };
        let s: String = [letter(xc, zc), letter(xt, zt)].iter().collect();
        // CNOT.
        let mut w: New = s.as_str().into();
        w.cnot(0, 1);
        let r = conj_cnot(P2 {
            phase: 0,
            xc,
            zc,
            xt,
            zt,
        });
        assert_eq!(
            (w.x_bit(0), w.z_bit(0), w.x_bit(1), w.z_bit(1)),
            (r.xc, r.zc, r.xt, r.zt),
            "CNOT bit map {s}"
        );
        // CZ.
        let mut w: New = s.as_str().into();
        w.cz(0, 1);
        let r = conj_cz(P2 {
            phase: 0,
            xc,
            zc,
            xt,
            zt,
        });
        assert_eq!(
            (w.x_bit(0), w.z_bit(0), w.x_bit(1), w.z_bit(1)),
            (r.xc, r.zc, r.xt, r.zt),
            "CZ bit map {s}"
        );
    }
}

// ===========================================================================
// Symplectic isometries (mirrors lean/PPVM/Pauli/Symplectic.lean)
// ===========================================================================

#[test]
fn clifford_generators_preserve_symplectic_form() {
    // hAct/sAct/cnotAct/czAct isometry: ω(Gv, Gw) = ω(v, w), run on the crate's
    // own Clifford map. Randomized n-qubit with random targets.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16] {
            for _ in 0..16 {
                let v0: New = random_pauli_string(&mut rng, n).as_str().into();
                let w0: New = random_pauli_string(&mut rng, n).as_str().into();
                let base = omega(&v0, &w0);
                let q = rng.random_range(0..n);
                let mut t = rng.random_range(0..n);
                while t == q {
                    t = rng.random_range(0..n);
                }

                for gate in 0..4 {
                    let (mut v, mut w) = (v0.clone(), w0.clone());
                    match gate {
                        0 => {
                            v.h(q);
                            w.h(q);
                        }
                        1 => {
                            v.s(q);
                            w.s(q);
                        }
                        2 => {
                            v.cnot(q, t);
                            w.cnot(q, t);
                        }
                        _ => {
                            v.cz(q, t);
                            w.cz(q, t);
                        }
                    }
                    assert_eq!(
                        omega(&v, &w),
                        base,
                        "seed {seed} gate {gate}: not an isometry"
                    );
                }
            }
        }
    }
}
