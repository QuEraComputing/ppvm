// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use fxhash::FxHashMap;
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::{HashFinalize, PauliStorage};
use std::f64::consts::PI;
use std::hash::BuildHasher;

use super::group::TranslationGroup;

impl TranslationGroup {
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
    assert_eq!(
        basis.len(),
        coeffs.len(),
        "basis and coeffs length mismatch"
    );
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
                    rep,
                    expected: *rep_coeff,
                    got_implied: implied_rep_coeff,
                    offending_pauli: *p,
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
