// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Basis truncation and enrichment: the magnitude prune, the rank cap,
//! and the capped leakage admission shared by every `pc_step` variant.
//!
//! All three are generic over the coefficient scalar ([`Coeff`]), so the
//! real adaptive path and the complex orbit-rep path run the same code.

use crate::scalar::Coeff;
use crate::word::Word;
use fxhash::{FxHashMap, FxHashSet};

/// Basis indices in descending coefficient magnitude. Leakage
/// accumulation walks the basis in this order so the running room-cap
/// keeps the entries most likely to be the true largest contributors.
pub(crate) fn order_by_desc_mag<T: Coeff>(coeffs: &[T]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..coeffs.len()).collect();
    order.sort_by(|&a, &b| desc_by_mag(coeffs[a], coeffs[b]));
    order
}

/// Keep only the `room` largest-magnitude entries of a live leakage
/// candidate map — `room` being the number of strings we could actually
/// admit to the basis, so there is no point tracking more. Applied after
/// each accumulation chunk.
pub(crate) fn cap_map_to_room<T: Coeff>(merged: &mut FxHashMap<Word, T>, room: usize) {
    if merged.len() <= room {
        return;
    }
    if room == 0 {
        merged.clear();
        return;
    }
    let mut mags: Vec<f64> = merged.values().map(|v| v.mag()).collect();
    let k = room.min(mags.len() - 1);
    let cutoff = nth_largest(&mut mags, k);
    merged.retain(|_, v| v.mag() >= cutoff);
}

/// Compact `basis` / `coeffs` in place: drop entries whose coefficient
/// magnitude is below `drop_tol` unless the word appears in `protected`.
/// No-op when `drop_tol ≤ 0`.
pub(crate) fn prune_basis<T: Coeff>(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<T>,
    drop_tol: f64,
    protected: &[Word],
) {
    if drop_tol <= 0.0 {
        return;
    }
    debug_assert_eq!(basis.len(), coeffs.len());
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    retain_in_place(basis, coeffs, |w, c| {
        c.mag() >= drop_tol || protected_set.contains(w)
    });
}

/// Global max-basis cap (PauliStrings.jl-style top-M trim): keep only the
/// `max_basis` largest-magnitude terms (protected strings always kept),
/// dropping the rest. Rank-based total-basis bound; dual of `drop_tol`.
/// A `max_basis` large enough to cover the whole basis is a no-op.
pub(crate) fn cap_basis<T: Coeff>(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<T>,
    max_basis: usize,
    protected: &[Word],
) {
    if basis.len() <= max_basis {
        return;
    }
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    let n_prot = basis.iter().filter(|w| protected_set.contains(w)).count();
    let slots = max_basis.saturating_sub(n_prot);
    let mut mags: Vec<f64> = basis
        .iter()
        .zip(coeffs.iter())
        .filter(|(w, _)| !protected_set.contains(w))
        .map(|(_, c)| c.mag())
        .collect();
    let cutoff = if slots == 0 {
        f64::INFINITY
    } else if slots >= mags.len() {
        return;
    } else {
        nth_largest(&mut mags, slots - 1)
    };
    retain_in_place(basis, coeffs, |w, c| {
        protected_set.contains(w) || c.mag() >= cutoff
    });
}

/// Add the largest leakage strings to the basis, up to the available room
/// `room = max_basis − basis.len()` — so the in-step basis (hence the
/// expm/leakage peak memory) never exceeds `max_basis`. New strings get
/// coefficient 0; the surrounding expm fills them. No magnitude filter: the
/// top-`room` by `|leakage|` are added (a large `max_basis` adds them all).
pub(crate) fn add_leakage_capped<T: Coeff>(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<T>,
    mut leak: Vec<(Word, T)>,
    max_basis: usize,
) {
    let room = max_basis.saturating_sub(basis.len());
    if leak.len() > room {
        if room > 0 {
            leak.select_nth_unstable_by(room - 1, |a, b| desc_by_mag(a.1, b.1));
        }
        leak.truncate(room);
    }
    for (w, _) in leak {
        basis.push(w);
        coeffs.push(T::zero());
    }
}

/// Keep the `basis`/`coeffs` entries satisfying `keep`, preserving order,
/// by swapping survivors down and truncating.
fn retain_in_place<T>(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<T>,
    mut keep: impl FnMut(&Word, &T) -> bool,
) {
    let mut write = 0;
    for read in 0..basis.len() {
        if keep(&basis[read], &coeffs[read]) {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}

/// The `k`-th largest element of `mags` (0-indexed), via a partial sort.
/// Reorders `mags`. Panics if `k >= mags.len()`.
fn nth_largest(mags: &mut [f64], k: usize) -> f64 {
    mags.select_nth_unstable_by(k, |a, b| {
        b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
    });
    mags[k]
}

/// Descending comparison by magnitude, NaN-tolerant.
fn desc_by_mag<T: Coeff>(a: T, b: T) -> std::cmp::Ordering {
    b.mag()
        .partial_cmp(&a.mag())
        .unwrap_or(std::cmp::Ordering::Equal)
}
