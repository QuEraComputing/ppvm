// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Kossakowski-form dissipator.
//!
//! For a family of operators `A_n` and a Hermitian, positive-semidefinite
//! pair matrix `K`, the adjoint dissipator is
//!
//! ```text
//! D*(O) = Σ_{n,m} K_nm ( A_n† O A_m − ½ {A_n† A_m, O} ).
//! ```
//!
//! This is the general GKSL form; the diagonal `K = diag(γ_k)` case is the
//! jump-operator form handled by [`crate::spec::JumpKind::General`].
//!
//! Each `(n, m)` pair is compiled once into a [`Pair`]. Hermitian-conjugate
//! pairs `(n,m)` and `(m,n)` are *folded* into a single upper-triangle entry:
//! both sandwiches produce the same output words with conjugate phase and the
//! final action keeps only the real part, so one entry that doubles-and-takes-
//! `Re` suffices, halving the pair count. [`PairShape`] records which of the
//! two a compiled pair is, and carries the extra term list that only the
//! folded case needs.

use crate::Error;
use crate::algebra::{
    COEFF_DROP_TOL, PauliTerm, comm_product, pauli_mul, phase_factor, precompute_adag_b,
    support_mask,
};
use crate::word::{Chunk, W_CHUNKS, Word, parse_pauli_string};
use fxhash::FxHashMap;
use num::Complex;
use std::collections::BTreeSet;

/// Relative tolerance for the Hermiticity check on `K`.
const HERMITICITY_TOL: f64 = 1e-10;

/// Sandwich table of a pair, grouped by the left word: one
/// `(P_a, [(P_b, coeff), …])` group per distinct `P_a`, so `P_a · p` is
/// computed once per group and reused across its `P_b` partners.
type SandwichGroups = Vec<(Word, Vec<(Word, Complex<f64>)>)>;

/// Which of the two compiled pair forms a [`Pair`] is.
///
/// The distinction changes the meaning of [`Pair::dd`] and selects the term
/// list used by the one-sided commutator path, so it is modelled as a sum
/// type rather than a flag: the off-diagonal-only term list cannot be
/// reached on a diagonal pair.
pub(crate) enum PairShape {
    /// `n == m`. [`Pair::dd`] is `K_nn · A_n†A_n`.
    Diagonal,
    /// `n < m`, folding in the conjugate `(m, n)` pair. [`Pair::dd`] is the
    /// Hermitian sum `2·Re(K_nm·A_n†A_m)` used by the both-sided
    /// anticommutator.
    OffDiagonal {
        /// The anti-Hermitian difference `−2i·Im(K_nm·A_n†A_m)`
        /// (pure-imaginary coefficients), used by the one-sided commutator
        /// of the folded conjugate pair.
        dd_anti: Vec<PauliTerm>,
    },
}

/// One compiled `(n, m)` pair of a Kossakowski dissipator.
pub(crate) struct Pair {
    sand: SandwichGroups,
    /// `A_n†A_m` scaled by `K_nm`; see [`PairShape`] for the exact form.
    dd: Vec<PauliTerm>,
    shape: PairShape,
    /// Support masks of `A_n` and `A_m`, for the one-sided fast path.
    left_mask: [Chunk; W_CHUNKS],
    right_mask: [Chunk; W_CHUNKS],
}

/// Compile a Kossakowski dissipator into one [`Pair`] per non-negligible
/// upper-triangle entry of `K`.
///
/// Returns the pairs alongside, for each pair, the union support of its two
/// operators, so the caller can index them by qubit.
pub(crate) fn compile(
    ops: &[Vec<(String, Complex<f64>)>],
    k: &[Vec<Complex<f64>>],
    n_qubits: usize,
) -> Result<Vec<(Pair, Vec<u32>)>, Error> {
    let max_abs = validate_k(k, ops.len())?;
    let (parsed, op_support) = parse_ops(ops, n_qubits)?;

    let pair_tol = COEFF_DROP_TOL * max_abs;
    let mut out = Vec::new();
    for n in 0..ops.len() {
        for m in n..ops.len() {
            if k[n][m].norm() <= pair_tol {
                continue;
            }
            let pair = compile_pair(&parsed[n], &parsed[m], k[n][m], n != m);
            let mut union: BTreeSet<u32> = op_support[n].iter().copied().collect();
            union.extend(op_support[m].iter().copied());
            out.push((pair, union.into_iter().collect()));
        }
    }
    Ok(out)
}

/// Check that `k` is square with side `n_ops` and Hermitian. Returns the
/// largest `|K_nm|`, which sets the scale for the negligible-pair cutoff.
fn validate_k(k: &[Vec<Complex<f64>>], n_ops: usize) -> Result<f64, Error> {
    if k.len() != n_ops {
        return Err(Error::LengthMismatch {
            what: "kossakowski ops and K rows",
            a: n_ops,
            b: k.len(),
        });
    }
    for (row, entries) in k.iter().enumerate() {
        if entries.len() != n_ops {
            return Err(Error::KMatrixRowLength {
                row,
                expected: n_ops,
                got: entries.len(),
            });
        }
    }

    let max_abs = k
        .iter()
        .flat_map(|row| row.iter().map(|c| c.norm()))
        .fold(0.0_f64, f64::max);

    // A non-Hermitian K is not a valid GKSL pair matrix and would produce an
    // action that does not preserve Hermiticity.
    let tol = HERMITICITY_TOL * max_abs.max(1.0);
    for (n, row_n) in k.iter().enumerate() {
        for (m, k_nm) in row_n.iter().enumerate().skip(n) {
            if (k_nm - k[m][n].conj()).norm() > tol {
                return Err(Error::KMatrixNotHermitian { n, m });
            }
        }
    }
    Ok(max_abs)
}

/// Parsed operator table: the Pauli terms of each `A_n`, and each `A_n`'s
/// union support.
type ParsedOps = (Vec<Vec<PauliTerm>>, Vec<Vec<u32>>);

/// Parse each operator's Pauli lincomb, returning the parsed terms and each
/// operator's union support.
fn parse_ops(ops: &[Vec<(String, Complex<f64>)>], n_qubits: usize) -> Result<ParsedOps, Error> {
    let mut parsed = Vec::with_capacity(ops.len());
    let mut supports = Vec::with_capacity(ops.len());
    for (i, op) in ops.iter().enumerate() {
        if op.is_empty() {
            return Err(Error::EmptyLincomb { index: i });
        }
        let mut terms = Vec::with_capacity(op.len());
        let mut union: BTreeSet<u32> = BTreeSet::new();
        for (s, c) in op {
            let (word, support) = parse_pauli_string(s, n_qubits)?;
            union.extend(support.iter().copied());
            terms.push(PauliTerm { word, coeff: *c });
        }
        parsed.push(terms);
        supports.push(union.into_iter().collect());
    }
    Ok((parsed, supports))
}

/// Compile the `(n, m)` entry with `A_n = a_terms`, `A_m = b_terms`.
fn compile_pair(
    a_terms: &[PauliTerm],
    b_terms: &[PauliTerm],
    k_nm: Complex<f64>,
    off_diag: bool,
) -> Pair {
    // A_n†A_m as `Σ γ_w W`, then scaled by K_nm.
    let adag_b = precompute_adag_b(a_terms, b_terms);
    let (dd, shape) = if off_diag {
        // Splitting K_nm·γ_w into its Hermitian and anti-Hermitian halves is
        // what lets the conjugate (m,n) pair be dropped: the (m,n) sandwich
        // contributes the complex conjugate, so the sum is 2·Re on the
        // both-sided path and 2i·Im on the one-sided one.
        let mut dd = Vec::with_capacity(adag_b.len());
        let mut dd_anti = Vec::with_capacity(adag_b.len());
        for t in &adag_b {
            let c = k_nm * t.coeff;
            if c.re.abs() > COEFF_DROP_TOL {
                dd.push(PauliTerm {
                    word: t.word,
                    coeff: Complex::new(2.0 * c.re, 0.0),
                });
            }
            if c.im.abs() > COEFF_DROP_TOL {
                dd_anti.push(PauliTerm {
                    word: t.word,
                    coeff: Complex::new(0.0, -2.0 * c.im),
                });
            }
        }
        (dd, PairShape::OffDiagonal { dd_anti })
    } else {
        let dd = adag_b
            .iter()
            .map(|t| PauliTerm {
                word: t.word,
                coeff: k_nm * t.coeff,
            })
            .collect();
        (dd, PairShape::Diagonal)
    };

    let sand = a_terms
        .iter()
        .map(|a| {
            let rights = b_terms
                .iter()
                .map(|b| (b.word, a.coeff.conj() * b.coeff * k_nm))
                .collect();
            (a.word, rights)
        })
        .collect();

    Pair {
        sand,
        dd,
        shape,
        left_mask: support_mask(a_terms),
        right_mask: support_mask(b_terms),
    }
}

impl Pair {
    /// Accumulate this pair's contribution to `L*(p)` into `local`.
    pub(crate) fn accumulate(&self, p: &Word, local: &mut FxHashMap<Word, Complex<f64>>) {
        let mut p_bits = [0 as Chunk; W_CHUNKS];
        for (i, slot) in p_bits.iter_mut().enumerate() {
            *slot = p.xbits.data[i] | p.zbits.data[i];
        }
        let hits = |mask: &[Chunk; W_CHUNKS]| (0..W_CHUNKS).any(|i| mask[i] & p_bits[i] != 0);
        let (hit_l, hit_r) = (hits(&self.left_mask), hits(&self.right_mask));

        // The pair is only visited when `p` overlaps at least one side, so
        // "not both" means exactly one.
        if hit_l && hit_r {
            self.accumulate_both_sided(p, local);
        } else {
            self.accumulate_one_sided(p, hit_r, local);
        }
    }

    /// One-sided fast path: when `p` is disjoint from one of the two
    /// operators the sandwich and anticommutator collapse to a commutator.
    /// For a diagonal pair with `D = K·A_n†A_m`:
    ///
    /// ```text
    /// p disjoint from A_n (left):  C = −½ [D, p]
    /// p disjoint from A_m (right): C = +½ [D, p]
    /// ```
    ///
    /// For a folded off-diagonal pair the two conjugate one-sided
    /// contributions combine into `±½ [F, p]` with the anti-Hermitian
    /// `F = dd_anti` and the *opposite* sign. With `[P_c, p] = −i·eps·out`
    /// from [`comm_product`], the term coefficient is `∓ t_c · (i/2) · eps`.
    fn accumulate_one_sided(
        &self,
        p: &Word,
        hit_r: bool,
        local: &mut FxHashMap<Word, Complex<f64>>,
    ) {
        let zero = Complex::new(0.0, 0.0);
        let (terms, half_i) = match &self.shape {
            PairShape::Diagonal => (
                &self.dd,
                if hit_r {
                    Complex::new(0.0, 0.5)
                } else {
                    Complex::new(0.0, -0.5)
                },
            ),
            PairShape::OffDiagonal { dd_anti } => (
                dd_anti,
                if hit_r {
                    Complex::new(0.0, -0.5)
                } else {
                    Complex::new(0.0, 0.5)
                },
            ),
        };
        for t in terms {
            let (out, eps) = comm_product(&t.word, p);
            if eps != 0.0 {
                *local.entry(out).or_insert(zero) += t.coeff * half_i * eps;
            }
        }
    }

    /// Both sides hit: full sandwich plus anticommutator. The sandwich is
    /// grouped by the left word so `P_a · p` is computed once per distinct
    /// `P_a` and reused across all its `P_b` partners. For a folded
    /// off-diagonal pair the sandwich is doubled and its real part taken
    /// (the conjugate `(m,n)` pair supplies the other half).
    fn accumulate_both_sided(&self, p: &Word, local: &mut FxHashMap<Word, Complex<f64>>) {
        let zero = Complex::new(0.0, 0.0);
        let fold = matches!(self.shape, PairShape::OffDiagonal { .. });
        for (wa, rights) in &self.sand {
            let (r_ap, phi1) = pauli_mul(wa, p);
            // Hoisted out of the inner loop: the fold is a property of the
            // pair, not of the term.
            if fold {
                for (wb, c0) in rights {
                    let (s, phi2) = pauli_mul(&r_ap, wb);
                    let v = c0 * phase_factor(phi1 + phi2);
                    *local.entry(s).or_insert(zero) += Complex::new(2.0 * v.re, 0.0);
                }
            } else {
                for (wb, c0) in rights {
                    let (s, phi2) = pauli_mul(&r_ap, wb);
                    *local.entry(s).or_insert(zero) += c0 * phase_factor(phi1 + phi2);
                }
            }
        }

        // −½{D, p}. For Pauli words, {P_c, p} = 2·sign·R when they commute
        // (P_c·p = sign·R) and 0 when they anti-commute; the ½ cancels the 2.
        for t in &self.dd {
            let (r, phase) = pauli_mul(&t.word, p);
            if phase & 1 == 0 {
                let sign = if phase == 0 { 1.0 } else { -1.0 };
                *local.entry(r).or_insert(zero) -= t.coeff * Complex::new(sign, 0.0);
            }
        }
    }
}
