// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::sum::PauliSum;
use fxhash::FxHashMap;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::Config;
use ppvm_traits::{HashFinalize, PauliStorage};
use std::hash::BuildHasher;

use super::group::TranslationGroup;

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
    assert_eq!(
        basis.len(),
        coeffs.len(),
        "basis and coeffs length mismatch"
    );
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
