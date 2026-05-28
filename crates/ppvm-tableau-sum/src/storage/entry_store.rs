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
    fn insert_or_merge_batch(
        &mut self,
        branches: Vec<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>,
        cutoff: &T::Coeff,
    ) -> bool;

    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&GeneralizedTableau<T, I, C>, &T::Coeff) -> bool;

    /// Rebuild identity caches after in-place tableau mutations and coalesce
    /// structurally equal entries by summing their probabilities. Returns true
    /// when at least one pair of entries was merged.
    fn merge_equal_entries(&mut self) -> bool;

    /// Reset `is_lost[addr0]` only on entries where it is currently set, update
    /// the corresponding loss-fingerprint delta, and coalesce any entries made
    /// structurally equal by that reset. Returns true when at least one pair of
    /// entries was merged.
    fn reset_loss_and_merge(&mut self, addr0: usize) -> bool;
}
