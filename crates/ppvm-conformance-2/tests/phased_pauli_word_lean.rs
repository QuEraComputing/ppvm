// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lean-oracle property tests for `ppvm-phased-pauli-word-2::PhasedPauliWord`:
//! the machine-checked semantics of `lean/PPVM/Pauli/**` reproduced as Rust
//! tests, run on the *actual* phased type (which — unlike the bare word — carries
//! the phase, so the sign half of every identity is exercised).
//!
//! Coverage (exhaustive for finite single-qubit cases, randomized for n-qubit):
//!
//! * `Phase.lean` — the single-qubit phased product **is** the group `𝒫₁`:
//!   `(φ,v)·(ψ,w) = (φ+ψ+phaseExp(v,w), v⊕w)`. We check the `Group` laws
//!   (`mul_assoc'` associativity, `one_mul'`/`mul_one'` identity,
//!   `inv_mul_cancel'` inverse) and `phaseExp_self` (`P·P = +I` up to the explicit
//!   sign) on the real `PhasedPauliWord`.
//! * `Word.lean` `phaseExpN_cocycle` / `phaseExpN_sub_comm` — the n-qubit twisted
//!   product is associative (2-cocycle) and `P·Q = (−1)^{ω} Q·P`.
//! * `Matrix.lean` `pauliMat_mul` / `tensorPauli_mul` — the accumulated phase is
//!   the base-`i` exponent of the genuine `ℤ[i]` matrix product of the two phased
//!   operators (single-qubit exhaustive over all phases; n-qubit randomized).
//! * `Conjugation.lean` `conjH_Y` (`HYH = −Y`), `conjS_X`/`conjS_Y`,
//!   `conjCNOT_sign`/`conjCZ_sign` — the Clifford conjugation **signs** the bare
//!   word dropped, now recovered by `Phased`'s real `PhaseTrack`, grounded in
//!   genuine `G·P·G†` matrices over `ℤ[i]` (single/two-qubit exhaustive; n-qubit
//!   randomized with random targets).
//! * the **extended** Clifford set (`S†`, `√X`, `√X†`, `√Y`, `√Y†`, `CY`) — no
//!   Lean `conj*` function exists for these, so they are grounded directly in
//!   `U† P U` over `ℤ[i]` with the square roots taken as the standard
//!   `exp(−i·π·P/4)` (scaled by `√2` to stay integral). This is the oracle both
//!   the `ppvm-traits-2` blanket (generator products) and `Phased`'s fused kernel
//!   are checked against.

use ppvm_conformance_2::{random_pauli_string, seeded_rng};
use ppvm_phased_pauli_word_2::PhasedPauliWord as Phased;
use ppvm_traits_2::{Clifford, PauliBits, Phase, Word};
use rand::RngExt;

const SEEDS: [u64; 12] = [1, 2, 3, 7, 11, 42, 99, 123, 777, 2024, 31337, 88888];
const PREFIXES: [&str; 4] = ["+", "+i", "-", "-i"];
const LETTERS: [char; 4] = ['I', 'X', 'Y', 'Z'];

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

/// The single-qubit Pauli `g(x,z) = iˣᶻ Xˣ Zᶻ` as a 2×2 ℤ[i] matrix, `Y = iXZ`.
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

/// Full n-fold Kronecker matrix of a bare Pauli string (qubit `0` leftmost).
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

/// The genuine matrix of a whole *phased* operator `i^φ · g(w)`.
fn phased_mat(w: &Phased) -> Mat {
    scalar_mul(ipow(w.phase().exponent()), &word_mat(&w.word().to_string()))
}

// ===========================================================================
// Small helpers on the phased type under test
// ===========================================================================

/// A `+`-phased identity word of width `n` (the group unit of `𝒫₁ⁿ`).
fn identity(n: usize) -> Phased {
    Phased::from(format!("+{}", "I".repeat(n)).as_str())
}

/// The group inverse `⟨-φ, x, z⟩` (each Pauli is an involution, so bits are kept
/// and only the phase negates; `Phase.lean inv` / `inv_mul_cancel'`).
fn group_inverse(w: &Phased) -> Phased {
    Phased::with_phase(*w.word(), w.phase().inverse())
}

/// Symplectic form `ω(p,q) = Σᵢ (x_pᵢ·z_qᵢ ⊕ z_pᵢ·x_qᵢ)` over 𝔽₂, read from the
/// phased words' inner bit accessors (matches `Symplectic.lean omegaFun`).
fn omega(p: &Phased, q: &Phased) -> u8 {
    let mut acc = 0u8;
    for i in 0..p.n_sites() {
        acc ^= (p.word().x_bit(i) && q.word().z_bit(i)) as u8;
        acc ^= (p.word().z_bit(i) && q.word().x_bit(i)) as u8;
    }
    acc
}

fn rand_phased(rng: &mut rand::rngs::StdRng, n: usize) -> Phased {
    let prefix = PREFIXES[rng.random_range(0..4usize)];
    Phased::from(format!("{prefix}{}", random_pauli_string(rng, n)).as_str())
}

// ===========================================================================
// 𝒫₁ group laws (Phase.lean) + n-qubit cocycle (Word.lean)
// ===========================================================================

#[test]
fn phased_product_is_group_mul_assoc() {
    // Phase.lean `mul_assoc'` / Word.lean `phaseExpN_cocycle`: the twisted product
    // is associative on both bits and folded phase.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let a = rand_phased(&mut rng, n);
            let b = rand_phased(&mut rng, n);
            let c = rand_phased(&mut rng, n);
            let left = &(&a * &b) * &c;
            let right = &a * &(&b * &c);
            assert_eq!(
                left.word().to_string(),
                right.word().to_string(),
                "assoc bits"
            );
            assert_eq!(left.phase(), right.phase(), "assoc phase seed {seed} n {n}");
        }
    }
}

#[test]
fn phased_product_identity_laws() {
    // Phase.lean `one_mul'` / `mul_one'`: the `+I` word is a two-sided identity.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let a = rand_phased(&mut rng, n);
            let id = identity(n);
            let left = &id * &a;
            let right = &a * &id;
            for (label, r) in [("one_mul", left), ("mul_one", right)] {
                assert_eq!(r.word().to_string(), a.word().to_string(), "{label} bits");
                assert_eq!(r.phase(), a.phase(), "{label} phase seed {seed}");
            }
        }
    }
}

#[test]
fn phased_product_inverse_law() {
    // Phase.lean `inv_mul_cancel'`: a · a⁻¹ = a⁻¹ · a = +I.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let a = rand_phased(&mut rng, n);
            let inv = group_inverse(&a);
            let id = identity(n);
            for (label, prod) in [("inv·a", &inv * &a), ("a·inv", &a * &inv)] {
                assert_eq!(
                    prod.word().to_string(),
                    id.word().to_string(),
                    "{label} bits"
                );
                assert_eq!(prod.phase(), Phase::Pos1, "{label} phase seed {seed}");
            }
        }
    }
}

#[test]
fn square_is_identity_up_to_sign() {
    // Phase.lean `phaseExp_self` (n-qubit `phaseExpN_self`): P·P = i^{2φ} · I, so
    // a `+`-phased word squares to exactly `+I`, and a general phase φ yields 2φ.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let a = rand_phased(&mut rng, n);
            let sq = &a * &a;
            let id = identity(n);
            assert_eq!(sq.word().to_string(), id.word().to_string(), "P² bits");
            // 2φ mod 4: the doubled explicit phase, residual is 0.
            let expected = Phase::from_exponent(2 * a.phase().exponent());
            assert_eq!(sq.phase(), expected, "P² phase seed {seed}");
        }
    }
}

#[test]
fn commutation_is_symplectic_form() {
    // Word.lean `phaseExpN_sub_comm`: (P·Q)·(Q·P)⁻¹ = (−1)^{ω(P,Q)} — bits agree,
    // phases differ by the symplectic form.
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &n in &[1usize, 2, 3, 5, 16, 60] {
            let p = rand_phased(&mut rng, n);
            let q = rand_phased(&mut rng, n);
            let pq = &p * &q;
            let qp = &q * &p;
            assert_eq!(
                pq.word().to_string(),
                qp.word().to_string(),
                "P·Q, Q·P bits"
            );
            let expected = if omega(&p, &q) == 0 {
                Phase::Pos1
            } else {
                Phase::Neg1
            };
            assert_eq!(
                pq.phase() * qp.phase().inverse(),
                expected,
                "seed {seed}: comm phase"
            );
        }
    }
}

// ===========================================================================
// Accumulated phase = base-i exponent of the ℤ[i] matrix product
// (Matrix.lean pauliMat_mul / tensorPauli_mul)
// ===========================================================================

#[test]
fn product_phase_equals_matrix_exponent_single_qubit_exhaustive() {
    // All 256 (phase,Pauli)×(phase,Pauli) single-qubit products vs the genuine
    // 2×2 ℤ[i] matrix product of the two phased operators.
    for pl in PREFIXES {
        for l in LETTERS {
            for pr in PREFIXES {
                for r in LETTERS {
                    let a = Phased::from(format!("{pl}{l}").as_str());
                    let b = Phased::from(format!("{pr}{r}").as_str());
                    let prod = &a * &b;
                    let lhs = matmul(&phased_mat(&a), &phased_mat(&b));
                    assert_eq!(lhs, phased_mat(&prod), "{pl}{l} · {pr}{r}");
                }
            }
        }
    }
}

#[test]
fn product_phase_equals_matrix_exponent_n_qubit_random() {
    // Same identity for n-fold Kronecker matrices (n ≤ 5 keeps 2ⁿ ≤ 32).
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 1..=5usize {
            let a = rand_phased(&mut rng, n);
            let b = rand_phased(&mut rng, n);
            let prod = &a * &b;
            let lhs = matmul(&phased_mat(&a), &phased_mat(&b));
            assert_eq!(lhs, phased_mat(&prod), "seed {seed}: {a} · {b}");
        }
    }
}

// ===========================================================================
// Clifford conjugation signs (Conjugation.lean) grounded in ℤ[i] matrices
// ===========================================================================

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

/// Conjugate transpose over `ℤ[i]` — needed because, unlike `H`/`CNOT`/`CZ`, the
/// extended Clifford generators are not Hermitian.
fn dagger(a: &Mat) -> Mat {
    let n = a.len();
    (0..n)
        .map(|i| (0..n).map(|j| a[j][i].conj()).collect())
        .collect()
}

// The extended Clifford gates over ℤ[i]. The square roots are defined as the
// standard `exp(-i·π·P/4)` (stim's `SQRT_X`/`SQRT_Y` up to global phase) and
// scaled by `√2` to stay integral, so `G† P G = 2·(U† P U)`.
fn sqrt2_sqrt_x() -> Mat {
    // √2·exp(−iπX/4) = I − iX.
    vec![vec![z(1, 0), z(0, -1)], vec![z(0, -1), z(1, 0)]]
}
fn sqrt2_sqrt_y() -> Mat {
    // √2·exp(−iπY/4) = I − iY.
    vec![vec![z(1, 0), z(-1, 0)], vec![z(1, 0), z(1, 0)]]
}
fn cy_gate() -> Mat {
    // |0⟩⟨0|⊗I + |1⟩⟨1|⊗Y, control = qubit 0.
    vec![
        vec![z(1, 0), z(0, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(1, 0), z(0, 0), z(0, 0)],
        vec![z(0, 0), z(0, 0), z(0, 0), z(0, -1)],
        vec![z(0, 0), z(0, 0), z(0, 1), z(0, 0)],
    ]
}

/// The single-qubit Pauli string for bit pair `(x, z)`.
fn letter(x: bool, zbit: bool) -> &'static str {
    match (x, zbit) {
        (false, false) => "I",
        (true, false) => "X",
        (false, true) => "Z",
        (true, true) => "Y",
    }
}

#[test]
fn single_qubit_conjugation_signs_grounded_in_zi_matrices() {
    // The phase `Phased`'s real PhaseTrack accumulates for H and S equals the
    // base-i exponent of the genuine matrix conjugation `G P G†` — closing the
    // "hand-derived sign" caveat.  H is Hermitian (H P H); S in this crate's
    // Heisenberg convention is `S† P S` (X ↦ −Y, the S/S† note in the crate's
    // `clifford.rs` and `Conjugation.lean`).
    for x in [false, true] {
        for zbit in [false, true] {
            let s = format!("+{}", letter(x, zbit));
            let mp = pauli_mat(x, zbit);

            // H: (√2 H) P (√2 H) = 2 H P H → expect 2·i^{delta}·result.
            let mut w = Phased::from(s.as_str());
            w.h(0);
            let lhs = matmul(&matmul(&sqrt2_h(), &mp), &sqrt2_h());
            assert_eq!(lhs, scalar_mul(z(2, 0), &phased_mat(&w)), "H sign {s}");

            // S: crate convention S† P S.
            let mut w = Phased::from(s.as_str());
            w.s(0);
            let lhs = matmul(&matmul(&s_dag(), &mp), &s_gate());
            assert_eq!(lhs, phased_mat(&w), "S sign {s}");
        }
    }
}

#[test]
fn two_qubit_conjugation_signs_grounded_in_zi_matrices() {
    // CNOT/CZ are Hermitian, self-inverse, integral → G P G, no scaling.  The
    // tracked phase matches the base-i exponent of the genuine matrix.
    for b in 0..16u8 {
        let (xc, zc, xt, zt) = (b & 1 != 0, b & 2 != 0, b & 4 != 0, b & 8 != 0);
        let s = format!("+{}{}", letter(xc, zc), letter(xt, zt));
        let mp = kron(&pauli_mat(xc, zc), &pauli_mat(xt, zt));

        let mut w = Phased::from(s.as_str());
        w.cnot(0, 1);
        let lhs = matmul(&matmul(&cnot_gate(), &mp), &cnot_gate());
        assert_eq!(lhs, phased_mat(&w), "CNOT sign {s}");

        let mut w = Phased::from(s.as_str());
        w.cz(0, 1);
        let lhs = matmul(&matmul(&cz_gate(), &mp), &cz_gate());
        assert_eq!(lhs, phased_mat(&w), "CZ sign {s}");
    }
}

/// The **extended** Clifford set (`S†`, `√X`, `√X†`, `√Y`, `√Y†`, `CY`) grounded
/// in genuine `ℤ[i]` matrix conjugation `U† P U`, the same way the generators
/// above are.
///
/// This is the oracle for the derivation in `ppvm-traits-2/src/pauli.rs`: the
/// blanket writes each extension gate as a *product* of audited generators
/// (`√X ≃ H·S·H`, …) and `Phased` supplies the fused equivalent, so the claim
/// that both realize the named gate is checked here against the gate's own
/// matrix — not against either implementation.
#[test]
fn extension_conjugation_signs_grounded_in_zi_matrices() {
    use ppvm_traits_2::CliffordExtensions;

    for x in [false, true] {
        for zbit in [false, true] {
            let s = format!("+{}", letter(x, zbit));
            let mp = pauli_mat(x, zbit);

            // S†: integral and unitary, so no scaling.
            let mut w = Phased::from(s.as_str());
            w.s_dag(0);
            let lhs = matmul(&matmul(&dagger(&s_dag()), &mp), &s_dag());
            assert_eq!(lhs, phased_mat(&w), "S† sign {s}");

            // The four square roots: G = √2·U, so G† P G = 2·(U† P U).
            for (name, g) in [
                ("√X", sqrt2_sqrt_x()),
                ("√X†", dagger(&sqrt2_sqrt_x())),
                ("√Y", sqrt2_sqrt_y()),
                ("√Y†", dagger(&sqrt2_sqrt_y())),
            ] {
                let mut w = Phased::from(s.as_str());
                match name {
                    "√X" => w.sqrt_x(0),
                    "√X†" => w.sqrt_x_dag(0),
                    "√Y" => w.sqrt_y(0),
                    _ => w.sqrt_y_dag(0),
                }
                let lhs = matmul(&matmul(&dagger(&g), &mp), &g);
                assert_eq!(lhs, scalar_mul(z(2, 0), &phased_mat(&w)), "{name} sign {s}");
            }
        }
    }

    // CY over the full two-qubit alphabet (integral and unitary → no scaling).
    for b in 0..16u8 {
        let (xc, zc, xt, zt) = (b & 1 != 0, b & 2 != 0, b & 4 != 0, b & 8 != 0);
        let s = format!("+{}{}", letter(xc, zc), letter(xt, zt));
        let mp = kron(&pauli_mat(xc, zc), &pauli_mat(xt, zt));

        let mut w = Phased::from(s.as_str());
        w.cy(0, 1);
        let lhs = matmul(&matmul(&dagger(&cy_gate()), &mp), &cy_gate());
        assert_eq!(lhs, phased_mat(&w), "CY sign {s}");
    }
}

#[test]
fn named_conjugation_sign_theorems() {
    // The specific `by decide` sign theorems of Conjugation.lean, asserted on the
    // real phased type.
    //
    // `conjH_Y`: HYH = −Y  (phase 2, bits unchanged).
    let mut w = Phased::from("+Y");
    w.h(0);
    assert_eq!(w.to_string(), "-Y", "conjH_Y");

    // `conjCNOT_sign`: the −1 delta appears exactly on XZ→−YY and YY→−XZ.
    let mut w = Phased::from("+XZ");
    w.cnot(0, 1);
    assert_eq!(w.to_string(), "-YY", "conjCNOT_sign XZ");
    let mut w = Phased::from("+YY");
    w.cnot(0, 1);
    assert_eq!(w.to_string(), "-XZ", "conjCNOT_sign YY");

    // `conjCZ_sign`: the −1 delta appears exactly on XY→−YX and YX→−XY.
    let mut w = Phased::from("+XY");
    w.cz(0, 1);
    assert_eq!(w.to_string(), "-YX", "conjCZ_sign XY");
    let mut w = Phased::from("+YX");
    w.cz(0, 1);
    assert_eq!(w.to_string(), "-XY", "conjCZ_sign YX");

    // `conjS_X`/`conjS_Y`: Conjugation.lean fixes S: X↦+Y, Y↦−X.  This crate's
    // Heisenberg S is the *adjoint* of that (the documented S/S† convention):
    // X↦−Y, Y↦+X.  We assert the crate's convention, which the matrix grounding
    // above confirms is a genuine `S† P S` conjugation.
    let mut w = Phased::from("+X");
    w.s(0);
    assert_eq!(w.to_string(), "-Y", "crate S: X ↦ −Y (adjoint of conjS_X)");
    let mut w = Phased::from("+Y");
    w.s(0);
    assert_eq!(w.to_string(), "+X", "crate S: Y ↦ +X (adjoint of conjS_Y)");
}

#[test]
fn n_qubit_conjugation_phase_grounded_in_zi_matrices() {
    // Randomized n-qubit (n ≤ 4 keeps 2ⁿ ≤ 16) with random targets: after a
    // single generator on a random phased word, the whole-operator matrix equals
    // `G M G†` (H √2-scaled per acting qubit).
    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for n in 2..=4usize {
            for _ in 0..8 {
                let w0 = rand_phased(&mut rng, n);
                let q = rng.random_range(0..n);
                let mut t = rng.random_range(0..n);
                while t == q {
                    t = rng.random_range(0..n);
                }
                let m0 = phased_mat(&w0);

                // H on q: conjugation by I⊗…⊗(√2H)_q⊗…, scaling the result by 2.
                let mut w = w0.clone();
                w.h(q);
                let g = embed_single(n, q, &sqrt2_h());
                let lhs = matmul(&matmul(&g, &m0), &g);
                assert_eq!(
                    lhs,
                    scalar_mul(z(2, 0), &phased_mat(&w)),
                    "seed {seed} H q{q}"
                );

                // S on q: S† M S (crate convention).
                let mut w = w0.clone();
                w.s(q);
                let gs = embed_single(n, q, &s_gate());
                let gsd = embed_single(n, q, &s_dag());
                let lhs = matmul(&matmul(&gsd, &m0), &gs);
                assert_eq!(lhs, phased_mat(&w), "seed {seed} S q{q}");
            }
        }
    }
}

/// Embed a single-qubit 2×2 gate `g` at qubit `q` in an `n`-qubit register as the
/// Kronecker `I⊗…⊗g⊗…⊗I` (qubit `0` leftmost).
fn embed_single(n: usize, q: usize, g: &Mat) -> Mat {
    let id = pauli_mat(false, false);
    let mut acc: Mat = vec![vec![z(1, 0)]];
    for i in 0..n {
        acc = kron(&acc, if i == q { g } else { &id });
    }
    acc
}
