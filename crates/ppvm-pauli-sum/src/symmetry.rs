// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Lattice translation symmetry groups for operator-space Pauli evolution.
//!
//! A [`TranslationGroup`] represents a finite abelian group `G` acting on
//! qubit positions by permutations. Given such a group, every Pauli word
//! belongs to a translation orbit, and operator dynamics that commute
//! with `G` can be tracked using **one canonical representative per
//! orbit** instead of all `|G|` orbit members — reducing per-step memory
//! and compute by a factor up to `|G|`.
//!
//! Following Teng, Chang, Rudolph, and Holmes (arXiv:2512.12094), this
//! module implements **plain (real-coefficient) merging** of Pauli sums
//! into orbit-representative form — see [`canonicalize_pauli_sum`] and
//! [`symmetry_merge_pauli_sum`]. This handles observables in the trivial
//! (`k=0`) symmetry sector, e.g. sums of single-Z operators over the
//! lattice.
//!
//! **Non-trivial momentum sectors (`k ≠ 0`)** are handled by
//! [`canonicalize_pauli_sum_complex`], which folds with the character
//! phase `χ_k(g)` of each translation. On the Python side, an operator in
//! sector `k` is carried as a *real pair* (real + imaginary components, two
//! real `PauliSum`s) and merged via `PauliSum.momentum_merge`, which reuses
//! this routine — letting gate-based Trotter evolution stay symmetry-
//! compressed in any momentum sector with real coefficients throughout.
//!
//! ## Data model
//!
//! A `TranslationGroup` is specified by a list of generator permutations
//! and their cyclic orders. The group order is the product of the orders.
//! For instance, a 2D `L × L` torus has two generators (translation in
//! x and y) each of order `L`.
//!
//! ## Canonicalization
//!
//! [`TranslationGroup::canonicalize`] returns the **lex-minimum** Pauli
//! word reachable from the input via group action. The ordering is the
//! standard `Ord` impl on `PauliWord` (compare `xbits`, then `zbits`).
//! All orbit members canonicalize to the same representative; orbits are
//! disjoint by construction, so the rep uniquely identifies the orbit.
//!
//! ## Merging
//!
//! [`canonicalize_pauli_sum`] takes parallel `Vec<Word>` / `Vec<f64>`
//! buffers (the representation used by ppvm-lindblad's adaptive
//! evolution) and replaces each Pauli by its canonical rep, summing
//! coefficients for collisions. The output is an orbit-rep basis with
//! coefficients equal to the sum of the input coefficients over each
//! orbit's members. For dynamics that commute with `G` and initial
//! states that are also `G`-invariant, this preserves the expectation
//! value of any `G`-invariant observable (paper's Theorem 1).
//!
//! See the dedicated tests for correctness against full-basis evolution
//! on small systems with no truncation.

use crate::sum::PauliSum;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::Config;
use ppvm_traits::{HashFinalize, PauliStorage, PauliWordTrait};
use fxhash::FxHashMap;
use num::Complex;
use std::f64::consts::PI;
use std::hash::BuildHasher;

/// A finite abelian symmetry group acting on qubit positions by
/// permutations.
///
/// Build via the convenience constructors [`Self::chain_1d`],
/// [`Self::torus_2d`], [`Self::torus_3d`], [`Self::ladder`], or
/// [`Self::from_generators`] for an arbitrary list of generator
/// permutations.
///
/// `perms[g]` is the permutation that **generator `g`** applies to qubit
/// indices: a qubit at position `q` moves to position `perms[g][q]`
/// under one application of generator `g`. `orders[g]` is the cyclic
/// order of generator `g` (i.e. applying it `orders[g]` times returns
/// the identity). The full group is the direct product of the cyclic
/// subgroups, with size `Π orders[g]`.
///
/// Only the **generators** are stored; the algorithm in
/// [`Self::canonicalize`] walks the group via mixed-radix increments.
#[derive(Debug, Clone)]
pub struct TranslationGroup {
    /// Number of qubits the group acts on.
    n_qubits: usize,
    /// One permutation per generator. `perms[g][q]` is the position
    /// that qubit `q` maps to under one application of generator `g`.
    perms: Vec<Vec<u32>>,
    /// Cyclic order of each generator.
    orders: Vec<u32>,
}

impl TranslationGroup {
    /// Construct from explicit generator permutations and orders.
    ///
    /// Each `perm` must be a permutation of `0..n_qubits`. Each `order`
    /// must satisfy `perm^order == identity`.
    pub fn from_generators(
        n_qubits: usize,
        perms: Vec<Vec<u32>>,
        orders: Vec<u32>,
    ) -> Self {
        assert_eq!(perms.len(), orders.len(), "perms and orders must match");
        for (g, perm) in perms.iter().enumerate() {
            assert_eq!(
                perm.len(),
                n_qubits,
                "generator {g} permutation has length {} != n_qubits {n_qubits}",
                perm.len()
            );
            let mut seen = vec![false; n_qubits];
            for &p in perm {
                assert!(
                    (p as usize) < n_qubits,
                    "generator {g} maps to out-of-range position {p}"
                );
                assert!(
                    !seen[p as usize],
                    "generator {g} is not a permutation (duplicate target {p})"
                );
                seen[p as usize] = true;
            }
        }
        Self {
            n_qubits,
            perms,
            orders,
        }
    }

    /// 1D chain of `n` sites with periodic boundary conditions.
    /// Single generator: cyclic shift by one site.
    pub fn chain_1d(n: usize) -> Self {
        let perm: Vec<u32> = (0..n).map(|q| ((q + 1) % n) as u32).collect();
        Self::from_generators(n, vec![perm], vec![n as u32])
    }

    /// 2D `lx × ly` torus, qubit at `(i, j)` indexed as `j*lx + i`.
    /// Two generators: x-shift (i → i+1 mod lx) and y-shift (j → j+1 mod ly).
    pub fn torus_2d(lx: usize, ly: usize) -> Self {
        let n = lx * ly;
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                (j * lx + (i + 1) % lx) as u32
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let (i, j) = (q % lx, q / lx);
                (((j + 1) % ly) * lx + i) as u32
            })
            .collect();
        Self::from_generators(n, vec![perm_x, perm_y], vec![lx as u32, ly as u32])
    }

    /// 3D `lx × ly × lz` torus, qubit at `(i, j, k)` indexed as
    /// `k*lx*ly + j*lx + i`.
    pub fn torus_3d(lx: usize, ly: usize, lz: usize) -> Self {
        let n = lx * ly * lz;
        let perm_x: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (k * lx * ly + j * lx + (i + 1) % lx) as u32
            })
            .collect();
        let perm_y: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (k * lx * ly + ((j + 1) % ly) * lx + i) as u32
            })
            .collect();
        let perm_z: Vec<u32> = (0..n)
            .map(|q| {
                let i = q % lx;
                let j = (q / lx) % ly;
                let k = q / (lx * ly);
                (((k + 1) % lz) * lx * ly + j * lx + i) as u32
            })
            .collect();
        Self::from_generators(
            n,
            vec![perm_x, perm_y, perm_z],
            vec![lx as u32, ly as u32, lz as u32],
        )
    }

    /// Multi-leg ladder: `l` sites along the chain × `n_legs` legs.
    /// Single generator: cyclic shift along the chain direction (all
    /// legs simultaneously). Qubit at `(leg, j)` indexed as
    /// `leg * l + j`. No translation along the leg axis (legs are
    /// distinguished).
    pub fn ladder(l: usize, n_legs: usize) -> Self {
        let n = l * n_legs;
        let perm: Vec<u32> = (0..n)
            .map(|q| {
                let leg = q / l;
                let j = q % l;
                (leg * l + (j + 1) % l) as u32
            })
            .collect();
        Self::from_generators(n, vec![perm], vec![l as u32])
    }

    /// Number of qubits the group acts on.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Number of generators (rank of the group as an abelian product).
    pub fn n_generators(&self) -> usize {
        self.perms.len()
    }

    /// Total group order: `Π orders[g]`.
    pub fn order(&self) -> usize {
        self.orders.iter().map(|&o| o as usize).product()
    }

    /// Permutation associated with the `g`-th generator (one application).
    pub fn generator_perm(&self, g: usize) -> &[u32] {
        &self.perms[g]
    }

    /// Cyclic order of the `g`-th generator.
    pub fn generator_order(&self, g: usize) -> u32 {
        self.orders[g]
    }

    /// Apply a single generator's permutation to a Pauli word, returning
    /// the resulting word.
    ///
    /// For each qubit `q` of the input, the corresponding `(xbit, zbit)`
    /// pair is placed at position `perm[q]` of the output.
    fn apply_generator<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
        g: usize,
    ) -> PauliWord<A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        let perm = &self.perms[g];
        let mut out: PauliWord<A, S, R> = PauliWord::new(self.n_qubits);
        for q in 0..self.n_qubits {
            let xb = w.get_xbit(q);
            let zb = w.get_zbit(q);
            if xb {
                out.set_xbit(perm[q] as usize, true);
            }
            if zb {
                out.set_zbit(perm[q] as usize, true);
            }
        }
        out.rehash();
        out
    }

    /// Lex-min canonical representative of `w`'s translation orbit
    /// under this group. Walks the full group via mixed-radix counters,
    /// keeping the smallest word seen.
    ///
    /// Total cost: `O(|G| × n_qubits)` per call.
    pub fn canonicalize<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> PauliWord<A, S, R>
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        debug_assert_eq!(
            w.n_qubits(),
            self.n_qubits,
            "word and group must agree on n_qubits"
        );
        if self.perms.is_empty() {
            return w.clone();
        }
        // Mixed-radix counter `(c[0], c[1], …)` ranges over
        // `0..orders[0] × 0..orders[1] × …`. We track the "current"
        // word obtained by applying generator `g` once each time
        // `c[g]` increments; rolling over `c[g]` means we apply
        // generator `g` exactly `orders[g]` times (= identity), so
        // `cur` returns to the orbit member that had `c[g..]` as its
        // tail and `0` in slots 0..g.
        //
        // The simplest correct implementation just enumerates: for each
        // group element index, build the corresponding word from scratch
        // by applying the right number of each generator.
        let mut best = w.clone();
        let order = self.order();
        let mut idx = 0usize;
        while idx < order {
            // Decode `idx` to mixed-radix counter `c`
            let mut rem = idx;
            let mut counters: Vec<u32> = Vec::with_capacity(self.perms.len());
            for &o in &self.orders {
                counters.push((rem as u32) % o);
                rem /= o as usize;
            }
            // Construct the group element's permutation by composing
            // `generator g` applied `c[g]` times, for each g.
            // We do this lazily by iterating over qubits.
            let mut cur = w.clone();
            for (g, &c) in counters.iter().enumerate() {
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            if cur < best {
                best = cur;
            }
            idx += 1;
        }
        best
    }

    /// Lex-min canonical representative `r` of `w` together with the
    /// **mixed-radix counter** `c = (c_0, c_1, …)` of the group element
    /// `g` such that `g·r = w`.
    ///
    /// In other words: if `r = self.canonicalize(w)`, this returns
    /// `(r, c)` where applying generator `i` exactly `c[i]` times in
    /// sequence to `r` produces `w`. The counter is used to compute
    /// momentum phases by the phase-aware merge routines.
    ///
    /// Same `O(|G| × n_qubits)` cost as `canonicalize`.
    pub fn canonicalize_with_shift<A, S, const R: bool>(
        &self,
        w: &PauliWord<A, S, R>,
    ) -> (PauliWord<A, S, R>, Vec<u32>)
    where
        A: PauliStorage,
        S: BuildHasher + Clone + Default + HashFinalize,
    {
        debug_assert_eq!(w.n_qubits(), self.n_qubits);
        if self.perms.is_empty() {
            return (w.clone(), Vec::new());
        }
        let mut best = w.clone();
        let mut best_counter: Vec<u32> = vec![0; self.perms.len()];
        let order = self.order();
        for idx in 0..order {
            // Decode `idx` to mixed-radix counter.
            let mut rem = idx;
            let mut counter: Vec<u32> = Vec::with_capacity(self.perms.len());
            for &o in &self.orders {
                counter.push((rem as u32) % o);
                rem /= o as usize;
            }
            // Build the candidate by applying generator `g` exactly
            // `counter[g]` times.
            let mut cur = w.clone();
            for (g, &c) in counter.iter().enumerate() {
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            if cur < best {
                best = cur;
                // We need the counter such that g·best = w. The loop
                // above computed cur = g·w with counter, so w = g^{-1}·cur.
                // For abelian cyclic groups, g^{-1} = g^{order-1}, i.e.
                // the counter `(orders[g] - counter[g]) mod orders[g]`.
                best_counter = counter
                    .iter()
                    .zip(self.orders.iter())
                    .map(|(&c, &o)| (o - c) % o)
                    .collect();
            }
        }
        (best, best_counter)
    }

    /// Momentum-sector character `χ_k(g) = exp(i Σ_g 2π · k[g] · counter[g] / orders[g])`
    /// where `k[g] ∈ ℤ` is the integer momentum mode along generator `g`
    /// (the corresponding wavenumber is `2π · k[g] / orders[g]`).
    ///
    /// `k.len()` must equal `self.n_generators()`. The character of the
    /// identity element (`counter = [0, …]`) is `1`. For the trivial
    /// (`k = [0, …]`) sector all characters are `1` — phase-aware merging
    /// reduces to plain merging.
    pub fn character(&self, k_modes: &[i32], counter: &[u32]) -> Complex<f64> {
        debug_assert_eq!(k_modes.len(), self.perms.len());
        debug_assert_eq!(counter.len(), self.perms.len());
        let mut phase = 0.0_f64;
        for ((&k, &c), &o) in k_modes.iter().zip(counter.iter()).zip(self.orders.iter()) {
            phase += 2.0 * PI * (k as f64) * (c as f64) / (o as f64);
        }
        Complex::from_polar(1.0, phase)
    }

    /// Iterate over all group elements applied to `w`. Yields `|G|`
    /// Pauli words (including `w` itself for the identity element).
    pub fn orbit<'a, A, S, const R: bool>(
        &'a self,
        w: &'a PauliWord<A, S, R>,
    ) -> impl Iterator<Item = PauliWord<A, S, R>> + 'a
    where
        A: PauliStorage + 'a,
        S: BuildHasher + Clone + Default + HashFinalize + 'a,
    {
        let order = self.order();
        (0..order).map(move |idx| {
            let mut rem = idx;
            let mut cur = w.clone();
            for (g, &o) in self.orders.iter().enumerate() {
                let c = (rem as u32) % o;
                rem /= o as usize;
                for _ in 0..c {
                    cur = self.apply_generator(&cur, g);
                }
            }
            cur
        })
    }
}

/// Replace `(basis, coeffs)` in-place with the orbit-representative
/// form: each Pauli word becomes its canonical rep, and coefficients
/// of words that collapse to the same rep are summed.
///
/// Output length ≤ input length. Entries whose summed coefficient
/// equals zero exactly are *not* removed — caller should run a final
/// `drop_tol` prune if desired.
///
/// For dynamics that commute with `group` and initial states that are
/// `group`-invariant (i.e. in the trivial momentum sector), this
/// preserves all `G`-invariant expectation values.
pub fn canonicalize_pauli_sum<A, S, const R: bool>(
    basis: &mut Vec<PauliWord<A, S, R>>,
    coeffs: &mut Vec<f64>,
    group: &TranslationGroup,
) where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(basis.len(), coeffs.len(), "basis and coeffs length mismatch");
    let mut merged: FxHashMap<PauliWord<A, S, R>, f64> =
        FxHashMap::with_capacity_and_hasher(basis.len(), Default::default());
    for (w, &c) in basis.iter().zip(coeffs.iter()) {
        let rep = group.canonicalize(w);
        *merged.entry(rep).or_insert(0.0) += c;
    }
    basis.clear();
    coeffs.clear();
    basis.reserve(merged.len());
    coeffs.reserve(merged.len());
    for (w, c) in merged {
        basis.push(w);
        coeffs.push(c);
    }
}

/// Replace `(basis, complex_coeffs)` in-place with the orbit-rep form
/// **projected onto momentum sector `k_modes`**.
///
/// Each Pauli `p` is replaced by its canonical rep `r`; the contribution
/// is `(1/|G|) · χ_k(g) · c_p` where `g` is the group element such that
/// `g · r = p` and `χ_k(g) = exp(2πi · Σ_g k_modes[g] · counter[g] / orders[g])`.
///
/// If the input was already a momentum-`k_modes` eigenstate (i.e. the
/// coefficients satisfy `c_{g·p} = χ_k(g)⁻¹ · c_p` for every orbit),
/// the output is the orbit-rep coefficients of that state unchanged.
/// Otherwise the merge discards the components in other sectors —
/// use [`check_momentum_sector`] beforehand to validate.
///
/// For the `k_modes = [0, 0, …]` (trivial) sector this reduces to plain
/// [`canonicalize_pauli_sum`] (real coefficients work, but on complex
/// input the result is complex with vanishing imaginary part).
pub fn canonicalize_pauli_sum_complex<A, S, const R: bool>(
    basis: &mut Vec<PauliWord<A, S, R>>,
    coeffs: &mut Vec<Complex<f64>>,
    group: &TranslationGroup,
    k_modes: &[i32],
) where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(basis.len(), coeffs.len(), "basis and coeffs length mismatch");
    assert_eq!(
        k_modes.len(),
        group.n_generators(),
        "k_modes length {} != number of generators {}",
        k_modes.len(),
        group.n_generators()
    );
    let inv_g: f64 = 1.0 / (group.order() as f64);
    let mut merged: FxHashMap<PauliWord<A, S, R>, Complex<f64>> =
        FxHashMap::with_capacity_and_hasher(basis.len(), Default::default());
    for (w, &c) in basis.iter().zip(coeffs.iter()) {
        let (rep, cnt) = group.canonicalize_with_shift(w);
        let chi = group.character(k_modes, &cnt);
        let contrib = inv_g * chi * c;
        *merged.entry(rep).or_insert(Complex::new(0.0, 0.0)) += contrib;
    }
    basis.clear();
    coeffs.clear();
    basis.reserve(merged.len());
    coeffs.reserve(merged.len());
    for (w, c) in merged {
        basis.push(w);
        coeffs.push(c);
    }
}

/// Verify that a `(basis, complex_coeffs)` Pauli sum lies entirely in
/// the momentum sector `k_modes` under `group`.
///
/// Concretely: for every orbit represented in the basis, all members
/// must satisfy `c_{g·r} = χ_k(g)⁻¹ · c_r` for some choice of orbit-rep
/// coefficient `c_r`.
///
/// Returns `Ok(())` on pass; `Err(SectorCheckError)` on fail with the
/// offending orbit-rep, expected coefficient, and actual coefficient.
///
/// Use this on a user-supplied initial state before feeding it to a
/// phase-aware merging pipeline — silently projecting a wrongly-typed
/// input throws away meaningful physics.
pub fn check_momentum_sector<A, S, const R: bool>(
    basis: &[PauliWord<A, S, R>],
    coeffs: &[Complex<f64>],
    group: &TranslationGroup,
    k_modes: &[i32],
    tol: f64,
) -> Result<(), SectorCheckError<A, S, R>>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    assert_eq!(basis.len(), coeffs.len());
    assert_eq!(k_modes.len(), group.n_generators());

    // Group entries by orbit rep, picking the first-seen member as
    // reference and checking later members against it.
    let mut reference: FxHashMap<PauliWord<A, S, R>, (Complex<f64>, Vec<u32>)> =
        FxHashMap::default();
    for (p, &c) in basis.iter().zip(coeffs.iter()) {
        let (rep, cnt) = group.canonicalize_with_shift(p);
        let chi = group.character(k_modes, &cnt);
        // expected c_p given the rep coefficient c_r:
        //   c_p = χ_k(g)⁻¹ · c_r,  where p = g·r
        // equivalently, c_r = χ_k(g) · c_p (a rearrangement).
        let implied_rep_coeff = chi * c;
        if let Some((rep_coeff, _ref_cnt)) = reference.get(&rep) {
            if (implied_rep_coeff - rep_coeff).norm() > tol * rep_coeff.norm().max(1.0) {
                return Err(SectorCheckError {
                    rep: rep.clone(),
                    expected: *rep_coeff,
                    got_implied: implied_rep_coeff,
                    offending_pauli: p.clone(),
                    offending_coeff: c,
                    shift: cnt.clone(),
                });
            }
        } else {
            reference.insert(rep, (implied_rep_coeff, cnt));
        }
    }
    Ok(())
}

/// Detail report for a failed [`check_momentum_sector`].
pub struct SectorCheckError<A: PauliStorage, S, const R: bool> {
    /// Canonical orbit representative for which the check failed.
    pub rep: PauliWord<A, S, R>,
    /// Coefficient that the *first* basis entry implied for `rep`.
    pub expected: Complex<f64>,
    /// Coefficient that `offending_pauli` implies for `rep` under the
    /// purported momentum sector.
    pub got_implied: Complex<f64>,
    /// The basis entry whose coefficient is inconsistent with the
    /// expected `rep` value.
    pub offending_pauli: PauliWord<A, S, R>,
    /// Original coefficient of `offending_pauli` in the input basis.
    pub offending_coeff: Complex<f64>,
    /// Counter encoding the group element `g` such that
    /// `g · rep == offending_pauli`.
    pub shift: Vec<u32>,
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Debug for SectorCheckError<A, S, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SectorCheckError {{ rep: <Word>, expected: {:?}, got_implied: {:?}, \
             offending: <Word>, offending_coeff: {:?}, shift: {:?} }}",
            self.expected, self.got_implied, self.offending_coeff, self.shift,
        )
    }
}

impl<A: PauliStorage, S, const R: bool> std::fmt::Display for SectorCheckError<A, S, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "input not in target momentum sector: orbit rep expected c={:?}, but \
             orbit member (shift {:?}, coeff {:?}) implies c={:?}",
            self.expected, self.shift, self.offending_coeff, self.got_implied,
        )
    }
}

/// Symmetry-merge a [`PauliSum`] in place: each Pauli word becomes its
/// canonical orbit representative, and entries collapsing to the same
/// rep accumulate coefficients.
///
/// This is the Trotter-mode counterpart to [`canonicalize_pauli_sum`]
/// (which operates on the `Vec<Word>, Vec<f64>` representation used by
/// `ppvm-lindblad`'s adaptive evolution). Same semantics: preserves all
/// `G`-invariant expectation values when the dynamics commutes with
/// `group` and the initial state is `group`-invariant.
///
/// Generic over the [`Config`] but constrained to PauliWord-backed
/// representations (i.e. not the loss-aware variant) since
/// canonicalization needs raw `(xbit, zbit)` access.
pub fn symmetry_merge_pauli_sum<T, A, S, const R: bool>(
    psum: &mut PauliSum<T>,
    group: &TranslationGroup,
) where
    T: Config<PauliWordType = PauliWord<A, S, R>>,
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    psum.map_add(|word, coeff| (group.canonicalize(word), coeff.clone()));
}

#[cfg(test)]
mod tests {
    use super::*;

    type W = PauliWord<[u8; 1], fxhash::FxBuildHasher, true>;

    fn word(s: &str) -> W {
        W::from(s)
    }

    #[test]
    fn chain_1d_canonicalizes_via_cyclic_shift() {
        let g = TranslationGroup::chain_1d(4);
        // All cyclic shifts of "IIXY" should canonicalize to the same rep.
        let candidates = ["IIXY", "IXYI", "XYII", "YIIX"];
        let canon: Vec<W> = candidates.iter().map(|s| g.canonicalize(&word(s))).collect();
        for c in &canon[1..] {
            assert_eq!(*c, canon[0], "all cyclic shifts must canonicalize to same rep");
        }
    }

    #[test]
    fn chain_1d_canonicalize_is_lex_min() {
        let g = TranslationGroup::chain_1d(4);
        let canon = g.canonicalize(&word("YIIX"));
        let orbit: Vec<W> = g.orbit(&word("YIIX")).collect();
        let min = orbit.iter().min().unwrap();
        assert_eq!(canon, *min);
    }

    #[test]
    fn orbit_has_correct_size_for_chain() {
        let g = TranslationGroup::chain_1d(4);
        // "XIII" has orbit of size 4 (full chain).
        let orbit: Vec<W> = g.orbit(&word("XIII")).collect();
        assert_eq!(orbit.len(), 4);
        // "XIXI" has orbit of size 2 (period-2 invariant); 4 elements
        // total in the orbit iterator, but only 2 unique.
        let orbit: Vec<W> = g.orbit(&word("XIXI")).collect();
        assert_eq!(orbit.len(), 4); // iterator yields |G|, including duplicates
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 2);
    }

    #[test]
    fn torus_2d_canonicalize() {
        // 3x2 torus, 6 qubits.
        let g = TranslationGroup::torus_2d(3, 2);
        assert_eq!(g.n_qubits(), 6);
        assert_eq!(g.order(), 6);
        // X at (0,0) — orbit is all 6 single-X positions.
        let w = word("XIIIII");
        let orbit: Vec<W> = g.orbit(&w).collect();
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 6);
        // All canonicalize to the same rep.
        let canon = g.canonicalize(&w);
        for u in &unique {
            assert_eq!(g.canonicalize(u), canon);
        }
    }

    #[test]
    fn ladder_canonicalize() {
        // 2-leg ladder, L=3 → 6 qubits, group order 3 (no swap of legs).
        let g = TranslationGroup::ladder(3, 2);
        assert_eq!(g.n_qubits(), 6);
        assert_eq!(g.order(), 3);
        // X on leg 0 site 0: orbit = {(0,0), (0,1), (0,2)}, NOT including leg 1 sites.
        let w = word("XIIIII"); // qubit 0 = X
        let orbit: Vec<W> = g.orbit(&w).collect();
        assert_eq!(orbit.len(), 3);
        let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
        assert_eq!(unique.len(), 3);
        // The orbit should be {qubit 0=X, qubit 1=X, qubit 2=X} — all leg 0.
        let expected: std::collections::HashSet<W> =
            ["XIIIII", "IXIIII", "IIXIII"].iter().map(|s| word(s)).collect();
        assert_eq!(unique, expected);
    }

    #[test]
    fn canonicalize_pauli_sum_merges_orbit_members() {
        let g = TranslationGroup::chain_1d(4);
        let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
        let mut coeffs: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
        // All four collapse to one rep with coeff 1+2+3+4 = 10.
        assert_eq!(basis.len(), 1);
        assert!((coeffs[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn canonicalize_pauli_sum_keeps_distinct_orbits() {
        let g = TranslationGroup::chain_1d(4);
        // Two distinct orbits: {XIII, ...} (size 4) and {ZIII, ...} (size 4).
        let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("ZIII"), word("IZII")];
        let mut coeffs: Vec<f64> = vec![1.0, 1.0, 2.0, 2.0];
        canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
        assert_eq!(basis.len(), 2);
        // Coefficients should be {2.0, 4.0} in some order.
        let mut cs = coeffs.clone();
        cs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((cs[0] - 2.0).abs() < 1e-12);
        assert!((cs[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn canonicalize_with_shift_round_trip() {
        // For each cyclic shift of "IIXY" by `a` positions, the shift
        // counter returned should reproduce the original word when
        // applied to the canonical rep.
        let g = TranslationGroup::chain_1d(4);
        for src in ["IIXY", "IXYI", "XYII", "YIIX"] {
            let w = word(src);
            let (rep, cnt) = g.canonicalize_with_shift(&w);
            // Apply gen 0 `cnt[0]` times to rep, should equal w.
            let mut cur = rep;
            for _ in 0..cnt[0] {
                cur = g.apply_generator(&cur, 0);
            }
            assert_eq!(cur, w, "shift {cnt:?} doesn't reproduce {src}");
        }
    }

    #[test]
    fn character_trivial_sector_is_one() {
        let g = TranslationGroup::chain_1d(4);
        // k=0 mode → character is always 1.
        for cnt in [vec![0u32], vec![1u32], vec![2u32], vec![3u32]] {
            let chi = g.character(&[0], &cnt);
            assert!((chi - Complex::new(1.0, 0.0)).norm() < 1e-12);
        }
    }

    #[test]
    fn character_obeys_unit_modulus() {
        let g = TranslationGroup::chain_1d(4);
        for k in 0..4 {
            for a in 0..4 {
                let chi = g.character(&[k], &[a as u32]);
                assert!(
                    (chi.norm() - 1.0).abs() < 1e-12,
                    "|χ_{k}(T^{a})| should be 1, got {}",
                    chi.norm()
                );
            }
        }
    }

    #[test]
    fn momentum_zero_complex_merge_matches_real_merge() {
        // k=0 sector: complex merge with all-real input should give
        // real-valued orbit-rep coefficients equal to the plain
        // canonicalize_pauli_sum result.
        let g = TranslationGroup::chain_1d(4);
        let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
        let real_coeffs = vec![1.0, 2.0, 3.0, 4.0];

        let mut basis_real = basis.clone();
        let mut coeffs_real = real_coeffs.clone();
        canonicalize_pauli_sum(&mut basis_real, &mut coeffs_real, &g);

        let mut basis_c = basis.clone();
        let mut coeffs_c: Vec<Complex<f64>> =
            real_coeffs.iter().map(|&v| Complex::new(v, 0.0)).collect();
        canonicalize_pauli_sum_complex(&mut basis_c, &mut coeffs_c, &g, &[0]);

        // Plain merge sums all coefficients onto the single orbit-rep:
        // 1+2+3+4 = 10. Complex merge does the same with a 1/|G|
        // prefactor, so we expect 10/4 = 2.5 on the rep.
        assert_eq!(basis_real.len(), 1);
        assert_eq!(basis_c.len(), 1);
        assert!((coeffs_real[0] - 10.0).abs() < 1e-12);
        assert!((coeffs_c[0].re - 2.5).abs() < 1e-12);
        assert!(coeffs_c[0].im.abs() < 1e-12);
    }

    #[test]
    fn momentum_eigenstate_check_passes() {
        // O = Σ_j e^{ikj} Z_j for k = 2π/4 (mode 1) is a momentum-k
        // eigenstate. check_momentum_sector should accept.
        let g = TranslationGroup::chain_1d(4);
        let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let k_mode: i32 = 1;
        // Sector condition: c_{T^a p} = e^{-2πi k a / N} c_p.
        // Picking c_{Z_0} = 1: c_{Z_a} = e^{-2πi · 1 · a / 4} = (-i)^a.
        let coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * (k_mode as f64) * (a as f64) / 4.0))
            .collect();
        let res = check_momentum_sector(&basis, &coeffs, &g, &[k_mode], 1e-10);
        assert!(res.is_ok(), "valid k-eigenstate failed sector check: {res:?}");
    }

    #[test]
    fn momentum_eigenstate_check_fails_for_wrong_sector() {
        // Same eigenstate as above, but check against the wrong momentum.
        let g = TranslationGroup::chain_1d(4);
        let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
            .collect();
        // Check against k=0 (constant) — should fail.
        let res = check_momentum_sector(&basis, &coeffs, &g, &[0], 1e-10);
        assert!(
            res.is_err(),
            "k=1 eigenstate wrongly passed as k=0 sector"
        );
    }

    #[test]
    fn momentum_eigenstate_round_trip_merge_preserves_rep_coeff() {
        // Merge a k=1 eigenstate; the orbit-rep coefficient should be
        // unchanged (= 1.0 for our chosen normalization, picking
        // c_{Z_0} = 1).
        let g = TranslationGroup::chain_1d(4);
        let mut basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
        let mut coeffs: Vec<Complex<f64>> = (0..4_i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
            .collect();
        canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &g, &[1]);
        assert_eq!(basis.len(), 1);
        // The canonical rep of single-Z orbit is Z_0 (lex-min of
        // {ZIII, IZII, IIZI, IIIZ} is IIIZ since 'I' < 'Z' lex-wise on
        // the (xbits, zbits) tuple; let's just check we got a single
        // entry with norm 1.
        assert!(
            (coeffs[0].norm() - 1.0).abs() < 1e-10,
            "expected |c_rep|=1, got {}",
            coeffs[0].norm()
        );
    }

    /// Trotter-mode end-to-end check that `PauliSum::symmetry_merge`
    /// matches plain Trotter evolution post-canonicalized.
    ///
    /// Setup: n=4 qubit chain, PBC, XY rotations on each bond. Initial
    /// operator `O(0) = Σ_j Z_j` is translation-invariant.
    ///
    /// **dt must be tiny.** First-order Trotter on a chain with PBC is
    /// only translation-equivariant up to `O(dt^2)` (gate-order
    /// commutator errors are NOT themselves T-symmetric). The
    /// "merge-after-each-step" trajectory and the "merge-at-end"
    /// trajectory therefore diverge by an amount proportional to that
    /// Trotter error. We test in the dt → 0 limit where the divergence
    /// is below FP noise.
    #[test]
    fn pauli_sum_symmetry_merge_matches_plain_trotter() {
        use crate::config::indexmap::ByteFxHashF64;
        use crate::prelude::*;

        type Cfg = ByteFxHashF64<1>;

        let n: usize = 4;
        // Tiny dt — Trotter per-step error scales as dt^2 and shows up
        // as a translation-non-equivariant correction; we want it below
        // FP noise at the tolerance we assert below (1e-7).
        let dt = 1e-5_f64;
        let n_steps = 2usize;
        let group = TranslationGroup::chain_1d(n);

        // Total-Z initial: O(0) = Σ_j Z_j (translation-invariant).
        let mut o_u: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        let mut o_m: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
        for j in 0..n {
            let mut s: Vec<char> = vec!['I'; n];
            s[j] = 'Z';
            let st: String = s.into_iter().collect();
            o_u += (st.as_str(), 1.0);
            o_m += (st.as_str(), 1.0);
        }
        assert_eq!(o_u.len(), n);
        assert_eq!(o_m.len(), n);

        // Apply XY Trotter steps to both copies. With merging, call
        // symmetry_merge_pauli_sum after each step.
        for _ in 0..n_steps {
            for j in 0..n {
                let nxt = (j + 1) % n;
                o_u.rxx(j, nxt, dt);
                o_u.ryy(j, nxt, dt);
                o_m.rxx(j, nxt, dt);
                o_m.ryy(j, nxt, dt);
            }
            symmetry_merge_pauli_sum(&mut o_m, &group);
        }

        // Canonicalize the un-merged result once at the end.
        symmetry_merge_pauli_sum(&mut o_u, &group);

        // Compare as (word → coeff) maps, FP tolerance.
        let u: FxHashMap<_, f64> = o_u.iter().map(|(w, c)| (w.clone(), *c)).collect();
        let m: FxHashMap<_, f64> = o_m.iter().map(|(w, c)| (w.clone(), *c)).collect();
        assert_eq!(
            u.len(),
            m.len(),
            "post-merge basis sizes differ: u={} vs m={}",
            u.len(),
            m.len()
        );
        let mut max_diff = 0.0_f64;
        for (w, &cu) in &u {
            let cm = *m.get(w).unwrap_or_else(|| {
                panic!("rep present in u but not in m: {:?}", w);
            });
            max_diff = max_diff.max((cu - cm).abs());
        }
        // At dt = 1e-5 over 2 steps, accumulated Trotter
        // commutator-induced T-eq error is ~2·dt^2·|H|^2 ≈ 1e-9; we
        // assert 1e-7 to leave safety margin.
        assert!(
            max_diff < 1e-7,
            "Trotter with-merging diverged from without-merging: max |Δc| = {max_diff:e}"
        );
    }
}
