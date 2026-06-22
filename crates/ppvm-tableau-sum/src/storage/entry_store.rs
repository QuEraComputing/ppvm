// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_tableau::{data::GeneralizedTableau, sparsevec::SparseVector};
use ppvm_traits::config::Config;

/// One branch produced by a noise channel: its tableau, coefficient, and the
/// cached `(word_fingerprint, phase_loss_hash)` pair, so a merge can recompute
/// the full fingerprint (`word_fp ^ phase_loss`) without re-hashing the tableau.
pub type Branch<T, I, C> = (GeneralizedTableau<T, I, C>, <T as Config>::Coeff, u64, u64);

pub trait EntryStore<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>>: Clone {
    fn with_capacity(cap: usize) -> Self;
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn iter<'a>(&'a self) -> impl Iterator<Item = (&'a GeneralizedTableau<T, I, C>, &'a T::Coeff)>
    where
        T: 'a,
        I: 'a,
        C: 'a;

    /// Mutate each entry in place. `FnMut` so the closure can capture e.g.
    /// `&mut Vec<branches>` from noise channels.
    fn for_each_mut<F>(&mut self, f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut T::Coeff);

    /// Like [`for_each_mut`](Self::for_each_mut), but first ensures cached
    /// fingerprints are current, then passes each entry's `word_fingerprint`
    /// and `phase_loss_hash` to the closure. Noise channels use this so a
    /// spawned branch can inherit its parent's word-hash (invariant under
    /// X/Y/Z and `is_lost`) and incrementally update the phase/loss hash,
    /// instead of re-hashing the tableau.
    fn for_each_mut_with_keys<F>(&mut self, f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut T::Coeff, u64, u64);

    /// O(1) flag that gates have touched every tableau and cached fingerprints
    /// are stale. Implementations recompute lazily on the next branching call.
    fn mark_keys_dirty(&mut self);

    /// Merge each branch into an existing entry whose tableau is structurally
    /// equal within threshold, or push it as new if its coefficient exceeds
    /// `cutoff`. Returns true if any incoming branch was dropped (caller
    /// renormalizes).
    ///
    /// Each branch carries its parent's word-fingerprint and its own
    /// phase/loss hash (third and fourth tuple fields), so the full
    /// fingerprint is `word_fp ^ phase_loss` — no re-hashing of the tableau.
    fn insert_or_merge_batch(&mut self, branches: Vec<Branch<T, I, C>>, cutoff: &T::Coeff) -> bool;

    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&GeneralizedTableau<T, I, C>, &T::Coeff) -> bool;

    /// Remove every entry whose tableau matches `pred` and return it together
    /// with its cached fingerprint components, so the caller can mutate the
    /// drained entries and re-insert them via [`insert_or_merge_batch`]. For
    /// storage backends that don't cache `word_fp` and `phase_loss` separately
    /// (e.g. map-keyed buckets), the split is `(fp, 0)`; either component may
    /// be XORed against the change delta — the merge sees their XOR.
    fn drain_where<F>(&mut self, pred: F) -> Vec<Branch<T, I, C>>
    where
        F: FnMut(&GeneralizedTableau<T, I, C>) -> bool;
}
