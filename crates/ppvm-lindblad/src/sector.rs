// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Momentum sectors of a translation group, and the phase-aware
//! canonicalization that drives orbit-representative evolution.
//!
//! On the orbit-rep path the state lives entirely in **orbit-rep form**
//! throughout: the basis contains only canonical translation-orbit
//! representatives and the coefficients are complex (one per rep). The
//! dynamics `L*` is computed with **phase-aware action** — for each
//! output Pauli `q`, we canonicalize `q` to its orbit rep `r_q` with
//! shift counter `cnt_q`, and accumulate
//! `χ_k(g_{cnt_q}) · v · c_r · |orbit_r| / |orbit_{r_q}|` (where `v` is
//! the matrix element of `L*` between input rep `r` and output `q`).
//! [`Sector::canonicalize_phase`] is that step.
//!
//! Coefficients are in the **averaged** convention: `c_r` is the plain
//! coefficient of the rep word, as produced by
//! `canonicalize_pauli_sum_complex`. The character-weighted action is
//! naturally the generator in the *summing* convention
//! `ĉ_r = |orbit_r| · c_r`, which is where the orbit-size ratio comes
//! from; it is 1 whenever both orbits are free.
//!
//! The orbit-rep basis is ~`|G|`× smaller than the full-basis
//! representation, throughout the entire evolution.
//!
//! ## Limitations
//!
//! - Callers are responsible for ensuring the input basis is in
//!   orbit-rep form (i.e. each entry is the canonical representative of
//!   its translation orbit). Use [`canonicalize_basis_to_rep`] if
//!   needed.
//! - A [`Sector`] is fixed for the duration of one
//!   [`LindbladSpec::pc_step_orbit_rep`](crate::LindbladSpec::pc_step_orbit_rep)
//!   call. To compute a full site-resolved profile, call it once per
//!   momentum mode and inverse-Fourier the results.

use crate::Word;
use num::Complex;
use ppvm_pauli_sum::symmetry::TranslationGroup;

/// A momentum sector of a translation group: the group `G` together with
/// one integer mode index per generator. The wavenumber along generator
/// `g` is `2π · k_modes[g] / group.generator_order(g)`; `k_modes = [0, …]`
/// is the trivial sector.
///
/// The two halves are meaningless apart — every phase-aware routine
/// needs both — so they travel as one value.
#[derive(Clone, Copy)]
pub struct Sector<'a> {
    group: &'a TranslationGroup,
    k_modes: &'a [i32],
}

impl<'a> Sector<'a> {
    pub fn new(group: &'a TranslationGroup, k_modes: &'a [i32]) -> Self {
        Self { group, k_modes }
    }

    /// Canonicalize `q` to its orbit representative `r_q` and return it
    /// alongside the character phase `χ_k(g_{cnt_q})` of the group
    /// element that maps `q` to `r_q`, and the number of **distinct**
    /// members of that orbit. The phase weights the matrix element of
    /// `L*` when it is accumulated onto `r_q`; the orbit size converts
    /// between the two coefficient conventions (see
    /// [`Self::orbit_size`]).
    ///
    /// `None` when `q`'s orbit cannot carry this sector (its stabilizer
    /// is incompatible with `k`): the coefficient of such a rep is
    /// identically zero, so the term is dropped.
    #[inline]
    pub fn canonicalize_phase(&self, q: &Word) -> Option<(Word, Complex<f64>, usize)> {
        let (rep, counter, orbit_size) = self.group.canonicalize_in_sector(q, self.k_modes)?;
        let phase = self.group.character(self.k_modes, &counter);
        Some((rep, phase, orbit_size))
    }

    /// Number of **distinct** members of `w`'s translation orbit, or
    /// `None` if the orbit cannot carry this sector.
    ///
    /// This is the factor between the two orbit-rep coefficient
    /// conventions: the *averaged* one, in which `c_r` is the plain
    /// coefficient of the rep word (what `canonicalize_pauli_sum_complex`
    /// and this crate's public orbit-rep API use), and the *summing* one
    /// `ĉ_r = |orbit_r| · c_r` (what `momentum_merge_pauli_sum_pair`
    /// uses). It is `|G|` only for free orbits.
    #[inline]
    pub fn orbit_size(&self, w: &Word) -> Option<usize> {
        self.group
            .canonicalize_in_sector(w, self.k_modes)
            .map(|(_, _, orbit_size)| orbit_size)
    }
}

/// Replace each entry of `basis` with its canonical orbit
/// representative under `group`. Pure rewrite; coefficients are
/// untouched. Useful to enforce the orbit-rep invariant before calling
/// [`LindbladSpec::pc_step_orbit_rep`](crate::LindbladSpec::pc_step_orbit_rep).
///
/// Does NOT deduplicate — if multiple input entries collapse to the
/// same rep, both are kept (caller should run a merge afterwards).
pub fn canonicalize_basis_to_rep(basis: &mut [Word], group: &TranslationGroup) {
    for w in basis.iter_mut() {
        *w = group.canonicalize(w);
    }
}
