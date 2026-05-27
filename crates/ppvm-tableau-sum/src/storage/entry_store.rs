use num::Complex;
use ppvm_runtime::config::Config;
use ppvm_tableau::{data::GeneralizedTableau, sparsevec::SparseVector};

pub trait EntryStore<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>>: Clone {
    fn with_capacity(cap: usize) -> Self;
    fn len(&self) -> usize;

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
    /// word-fingerprints are current, then passes each entry's
    /// word-fingerprint to the closure. Noise channels use this so a spawned
    /// branch can inherit its parent's word-hash (invariant under X/Y/Z and
    /// `is_lost`) instead of re-hashing the Pauli words.
    fn for_each_mut_with_word_key<F>(&mut self, f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut T::Coeff, u64);

    /// O(1) flag that gates have touched every tableau and cached fingerprints
    /// are stale. Implementations recompute lazily on the next branching call.
    fn mark_keys_dirty(&mut self);

    /// Merge each branch into an existing entry whose tableau is structurally
    /// equal within threshold, or push it as new if its coefficient exceeds
    /// `cutoff`. Returns true if any incoming branch was dropped (caller
    /// renormalizes).
    ///
    /// Each branch carries its parent's word-fingerprint (third tuple field);
    /// since the branch ops leave `word` unchanged, the full fingerprint is
    /// recovered as `word_fp ^ phase_lost_fingerprint(branch)` without
    /// re-hashing the Pauli words.
    fn insert_or_merge_batch(
        &mut self,
        branches: Vec<(GeneralizedTableau<T, I, C>, T::Coeff, u64)>,
        cutoff: &T::Coeff,
    ) -> bool;

    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&GeneralizedTableau<T, I, C>, &T::Coeff) -> bool;
}
