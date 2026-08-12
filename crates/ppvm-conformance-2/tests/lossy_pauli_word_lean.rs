// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for `ppvm-lossy-pauli-word-2::LossyPauliWord`.
//!
//! Two kinds of facts (design: task brief; `word-data-structures.md` §"Lossy
//! Pauli word"; Lean spec `lean/PPVM/Pauli/**`):
//!
//! 1. **Pauli algebra on PRESENT (non-lost) sites** — the present-site projection
//!    of a lossy word (`Lost ↦ I`) is an ordinary Pauli word, so it must satisfy
//!    the very same machine-checked `PauliWord` oracles: `phaseExp` equals the
//!    ℤ[i] (Gaussian-integer) matrix exponent (`Matrix.lean pauliMat_mul`), the
//!    summed phase is a 2-cocycle (`Word.lean phaseExpN_cocycle`), `P·P = +I`
//!    (`phaseExpN_self`), `P·Q = (−1)^{ω} Q·P` (`phaseExpN_sub_comm`), and each
//!    Clifford generator is a symplectic isometry (`Symplectic.lean *Act_isometry`,
//!    run on the lossy word's *own* `Clifford`). These reuse the same ℤ[i] matrix
//!    reference `pauli_word_lean.rs` grounds `PauliWord` in.
//!
//! 2. **Loss-specific data-structure invariants** — genuinely new facts with no
//!    algebraic Lean counterpart, asserted as property tests:
//!    * a lost site has canonical bits `(x, z, lost) = (0, 0, 1)`;
//!    * loss is exclusive with `I/X/Y/Z` (a present site is never lost);
//!    * `loss_weight()` counts exactly the lost sites; and
//!    * loss is untouched by the Clifford bit ops (the lost-qubit set is
//!      invariant under H/S/CNOT/CZ).
//!
//! Exhaustive for the finite single-qubit cases, randomized for n-qubit.

use ppvm_conformance_2::{
    GateOp, LOSSY_LETTERS, PAULIS, random_circuit, random_lossy_pauli_string, seeded_rng,
};
use ppvm_lossy_pauli_word_2::LossyPauliWord as NewLossyWord;
use ppvm_pauli_word_2::PauliWord as NewWord;
use ppvm_traits_2::{Clifford, KeyProduct, LossySite, Pauli, PauliBits, Phase, Word};
use rand::RngExt;
use rand::rngs::StdRng;

type NewLossy = NewLossyWord<u64>;
type New = NewWord<u64>;

const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];

// ===========================================================================
// ℤ[i] matrix reference (mirrors lean/PPVM/Pauli/Matrix.lean; identical to the
// reference `pauli_word_lean.rs` grounds the ordinary word in)
// ===========================================================================

/// A Gaussian integer `ℤ[i]` (exact; no floats).
type Z = num::Complex<i64>;
/// A dense matrix over `ℤ[i]`.
type Mat = Vec<Vec<Z>>;

#[inline]
fn z(re: i64, im: i64) -> Z {
    num::Complex::new(re, im)
}

/// `iᵏ ∈ ℤ[i]` for exponent `k`.
fn ipow(k: u8) -> Z {
    match k & 3 {
        0 => z(1, 0),
        1 => z(0, 1),
        2 => z(-1, 0),
        _ => z(0, -1),
    }
}

/// The single-qubit Pauli `g(x,z) = iˣᶻ Xˣ Zᶻ` as a 2×2 ℤ[i] matrix (`Y = iXZ`).
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
// Present-site projection: Lost ↦ I; present ↦ its Pauli.
// ===========================================================================

/// The present-site projection of a lossy word as a plain Pauli string
/// (`Lost ↦ I`), which the ordinary-word oracles then act on.
fn present_projection(w: &NewLossy) -> String {
    (0..w.n_sites())
        .map(|i| match w.get(i) {
            LossySite::Present(Pauli::I) | LossySite::Lost => 'I',
            LossySite::Present(Pauli::X) => 'X',
            LossySite::Present(Pauli::Y) => 'Y',
            LossySite::Present(Pauli::Z) => 'Z',
        })
        .collect()
}

/// The ordinary `PauliWord` on the present-site projection of a lossy word.
fn projected_word(w: &NewLossy) -> New {
    present_projection(w).as_str().into()
}

// ===========================================================================
// phaseExp on the present projection == ℤ[i] matrix exponent
// ===========================================================================

const LETTERS: [char; 4] = ['I', 'X', 'Y', 'Z'];

#[test]
fn present_single_qubit_phase_matches_matrix_exhaustive() {
    // Every present single-qubit case (Lost projects to I and is covered by 'I').
    for p in LETTERS {
        for q in LETTERS {
            let a: New = p.to_string().as_str().into();
            let b: New = q.to_string().as_str().into();
            let (r, phase): (New, Phase) = a.key_mul(&b);
            let lhs = matmul(&pauli_mat_of_char(p), &pauli_mat_of_char(q));
            let rhs = scalar_mul(ipow(phase.exponent()), &word_mat(&r.to_string()));
            assert_eq!(lhs, rhs, "{p}·{q}: present-site phase != ℤ[i] matrix");
        }
    }
}

#[test]
fn present_projection_phase_matches_matrix_n_qubit_random() {
    // n ≤ 5 keeps 2ⁿ ≤ 32; project two random lossy words to their present sites
    // and check the twisted-product phase against the ℤ[i] matrix product.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=5usize {
            let lv: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            let lw: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            let pv = present_projection(&lv);
            let pw = present_projection(&lw);
            let a: New = pv.as_str().into();
            let b: New = pw.as_str().into();
            let (r, phase): (New, Phase) = a.key_mul(&b);

            let lhs = matmul(&word_mat(&pv), &word_mat(&pw));
            let rhs = scalar_mul(ipow(phase.exponent()), &word_mat(&r.to_string()));
            assert_eq!(lhs, rhs, "seed {seed}: present projection {pv}·{pw}");
        }
    }
}

// ===========================================================================
// cocycle / self / commutation on the present projection
// ===========================================================================

#[test]
fn present_projection_cocycle_associativity_random() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let u = projected_word(&random_lossy_pauli_string(&mut rng, n).as_str().into());
            let v = projected_word(&random_lossy_pauli_string(&mut rng, n).as_str().into());
            let w = projected_word(&random_lossy_pauli_string(&mut rng, n).as_str().into());

            let (uv, p1) = u.key_mul(&v);
            let (uv_w, p2) = uv.key_mul(&w);
            let left = p1.compose(p2);

            let (vw, p3) = v.key_mul(&w);
            let (u_vw, p4) = u.key_mul(&vw);
            let right = p3.compose(p4);

            assert_eq!(uv_w, u_vw, "seed {seed} n {n}: associativity bits");
            assert_eq!(left, right, "seed {seed} n {n}: cocycle phase");
        }
    }
}

#[test]
fn present_projection_square_is_plus_identity_random() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let l: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            let p = projected_word(&l);
            let (r, phase) = p.key_mul(&p);
            assert_eq!(r, New::new(n), "seed {seed}: P² word != I");
            assert_eq!(phase, Phase::Pos1, "seed {seed}: P² phase != +1");
        }
    }
}

/// Symplectic form `ω(p,q) = Σᵢ (x_pᵢ·z_qᵢ ⊕ z_pᵢ·x_qᵢ)` over 𝔽₂, read from the
/// **lossy** word's own bit accessors. A lost site has `x = z = 0`, so it never
/// contributes: `ω` over the lossy word equals `ω` over its present projection.
fn omega_lossy(p: &NewLossy, q: &NewLossy) -> u8 {
    let mut acc = 0u8;
    for i in 0..p.n_sites() {
        acc ^= (p.x_bit(i) && q.z_bit(i)) as u8;
        acc ^= (p.z_bit(i) && q.x_bit(i)) as u8;
    }
    acc
}

#[test]
fn present_projection_commutation_is_symplectic_form_random() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let lp: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            let lq: NewLossy = random_lossy_pauli_string(&mut rng, n).as_str().into();
            let p = projected_word(&lp);
            let q = projected_word(&lq);
            let (rpq, pq) = p.key_mul(&q);
            let (rqp, qp) = q.key_mul(&p);
            assert_eq!(rpq, rqp, "seed {seed}: P·Q, Q·P differ on bits");

            // ω read off the lossy words equals ω off the projections (lost = 0).
            let expected = if omega_lossy(&lp, &lq) == 0 {
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
// Clifford generators are symplectic isometries — run on the lossy Clifford
// ===========================================================================

/// A random lossy string with the given shared loss mask: `L` where `lost[i]`,
/// otherwise a random present Pauli. The shared mask is what makes the
/// loss-guarded Clifford act as the *same* symplectic map on two operands (the
/// guard's decision depends only on the loss plane), which is the precondition
/// for the isometry ω(Gv, Gw) = ω(v, w) to hold.
fn lossy_with_mask(rng: &mut StdRng, lost: &[bool]) -> String {
    lost.iter()
        .map(|&l| {
            if l {
                'L'
            } else {
                PAULIS[rng.random_range(0..4usize)]
            }
        })
        .collect()
}

#[test]
fn lossy_clifford_generators_preserve_symplectic_form() {
    // hAct/sAct/cnotAct/czAct isometry: ω(Gv, Gw) = ω(v, w), run on the lossy
    // word's *own* Clifford (which no-ops on lost qubits). Because lost sites
    // carry x = z = 0 they never contribute to ω, and a *shared* loss mask makes
    // the guarded gate the same Sp map on both operands, so the isometry holds
    // exactly as for an ordinary word on the present-site subspace.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16] {
            for _ in 0..16 {
                let lost: Vec<bool> = (0..n).map(|_| rng.random_range(0..5usize) == 4).collect();
                let v0: NewLossy = lossy_with_mask(&mut rng, &lost).as_str().into();
                let w0: NewLossy = lossy_with_mask(&mut rng, &lost).as_str().into();
                let base = omega_lossy(&v0, &w0);
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
                        omega_lossy(&v, &w),
                        base,
                        "seed {seed} gate {gate}: lossy Clifford not an isometry"
                    );
                }
            }
        }
    }
}

// ===========================================================================
// Loss-specific data-structure invariants
// ===========================================================================

#[test]
fn lost_site_has_canonical_bits_and_loss_is_exclusive() {
    // Exhaustive single-qubit alphabet, then randomized n-qubit.
    for c in LOSSY_LETTERS {
        let w: NewLossy = c.to_string().as_str().into();
        assert_site_invariants(&w, &c.to_string());
    }
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let s = random_lossy_pauli_string(&mut rng, n);
            let w: NewLossy = s.as_str().into();
            assert_site_invariants(&w, &s);
        }
    }
}

/// A lost site is `(x, z, lost) = (0, 0, 1)`; a present site is never lost and its
/// `(x, z)` bits encode its Pauli.
fn assert_site_invariants(w: &NewLossy, s: &str) {
    for i in 0..w.n_sites() {
        match w.get(i) {
            LossySite::Lost => {
                assert!(NewLossyWord::is_lost(w, i), "site {i} of {s} reports lost");
                assert!(!w.x_bit(i), "lost site {i} of {s} has x=0");
                assert!(!w.z_bit(i), "lost site {i} of {s} has z=0");
            }
            LossySite::Present(p) => {
                // Loss is exclusive with a present Pauli.
                assert!(
                    !NewLossyWord::is_lost(w, i),
                    "present site {i} of {s} is not lost"
                );
                let (x, zb) = match p {
                    Pauli::I => (false, false),
                    Pauli::X => (true, false),
                    Pauli::Z => (false, true),
                    Pauli::Y => (true, true),
                };
                assert_eq!((w.x_bit(i), w.z_bit(i)), (x, zb), "present bits {i} of {s}");
            }
        }
    }
}

#[test]
fn loss_weight_counts_lost_sites() {
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let s = random_lossy_pauli_string(&mut rng, n);
            let w: NewLossy = s.as_str().into();
            let by_symbol = s.chars().filter(|&c| c == 'L').count();
            let by_read = (0..n).filter(|&i| NewLossyWord::is_lost(&w, i)).count();
            assert_eq!(w.loss_weight(), by_symbol, "loss_weight vs symbols {s}");
            assert_eq!(w.loss_weight(), by_read, "loss_weight vs is_lost reads {s}");
        }
    }
}

#[test]
fn clifford_leaves_loss_mask_invariant() {
    // The lost-qubit set is invariant under H/S/CNOT (the bit ops no-op on loss).
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[2usize, 3, 5, 16] {
            let s = random_lossy_pauli_string(&mut rng, n);
            let mut w: NewLossy = s.as_str().into();
            let before: Vec<bool> = (0..n).map(|i| NewLossyWord::is_lost(&w, i)).collect();

            let circuit = random_circuit(&mut rng, n, 200);
            for op in circuit {
                match op {
                    GateOp::H(q) => w.h(q),
                    GateOp::S(q) => w.s(q),
                    GateOp::Cnot(c, t) => w.cnot(c, t),
                    GateOp::Rx(..) | GateOp::Rz(..) => {}
                }
            }

            let after: Vec<bool> = (0..n).map(|i| NewLossyWord::is_lost(&w, i)).collect();
            assert_eq!(
                before, after,
                "loss mask changed under Clifford {s} seed {seed}"
            );
        }
    }
}
