use std::ops::AddAssign;

use fxhash::{FxHashMap, FxHashSet};
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use smallvec::SmallVec;

use crate::storage::{
    EntryStore, loss_mask, phase_loss_hash, structurally_equal, word_fingerprint,
};

#[derive(Clone)]
pub struct VecStorage<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> {
    pub entries: Vec<(GeneralizedTableau<T, I, C>, T::Coeff)>,
    pub fingerprints: Vec<u64>,
    /// Cached `word_fingerprint` per entry, kept in lockstep with `entries`.
    /// A branch inherits its parent's value (its Pauli words are unchanged),
    /// so the merge avoids re-hashing the words — the dominant fingerprint cost.
    pub word_fingerprints: Vec<u64>,
    /// Cached `phase_loss_hash` per entry, kept in lockstep with `entries`.
    /// XOR-combinable, so a branch inherits its parent's value and updates only
    /// the rows it changed — the merge never walks the phases again.
    pub phase_loss_hashes: Vec<u64>,
    pub dirty: bool,
    /// Reusable scratch buffer for `structurally_equal`'s coefficient lookup
    /// map. Cleared and refilled per call; keeps its capacity across calls.
    scratch: FxHashMap<I, Complex<T::Coeff>>,
}

impl<T, I, C> VecStorage<T, I, C>
where
    T: Config,
    T::Coeff: One + Zero + Clone + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    fn insert_or_merge_entry(
        &mut self,
        tab: GeneralizedTableau<T, I, C>,
        p: T::Coeff,
        word_fp: u64,
        phase_loss: u64,
        fp_index: &mut FxHashMap<u64, Vec<usize>>,
        sum_cutoff: &T::Coeff,
    ) -> bool {
        // The branch carries both components — its inherited word-fingerprint
        // and its incrementally-maintained phase/loss hash — so the full
        // fingerprint is their XOR, with no walk over the tableau.
        let fp = word_fp ^ phase_loss;

        // Only run the full equality check on entries whose fingerprint matches.
        let mut found: Option<usize> = None;
        if let Some(candidates) = fp_index.get(&fp) {
            for &i in candidates {
                if structurally_equal(&self.entries[i].0, &tab, &mut self.scratch) {
                    found = Some(i);
                    break;
                }
            }
        }

        let mut needs_normalize = false;
        match found {
            Some(i) => {
                let p0 = &self.entries[i].1;
                self.entries[i].1 = p0.clone() + p;
            }

            None => {
                if &p > sum_cutoff {
                    let new_idx = self.entries.len();
                    self.entries.push((tab, p));
                    self.fingerprints.push(fp);
                    self.word_fingerprints.push(word_fp);
                    self.phase_loss_hashes.push(phase_loss);
                    fp_index.entry(fp).or_default().push(new_idx);
                } else {
                    needs_normalize = true;
                }
            }
        }

        needs_normalize
    }

    /// Recompute `fingerprints` and `word_fingerprints` from scratch when a
    /// Clifford/T gate has invalidated them. No-op when clean.
    fn rebuild_fingerprints_if_dirty(&mut self) {
        if !self.dirty {
            return;
        }
        self.fingerprints.clear();
        self.word_fingerprints.clear();
        self.phase_loss_hashes.clear();
        for (t, _) in self.entries.iter() {
            let wfp = word_fingerprint(t);
            let plh = phase_loss_hash(t);
            self.word_fingerprints.push(wfp);
            self.phase_loss_hashes.push(plh);
            self.fingerprints.push(wfp ^ plh);
        }
        self.dirty = false;
    }
}

impl<T, I, C> EntryStore<T, I, C> for VecStorage<T, I, C>
where
    T: Config,
    T::Coeff: One + Zero + Clone + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    fn with_capacity(cap: usize) -> Self {
        Self {
            entries: Vec::with_capacity(cap),
            fingerprints: Vec::with_capacity(cap),
            word_fingerprints: Vec::with_capacity(cap),
            phase_loss_hashes: Vec::with_capacity(cap),
            dirty: false,
            scratch: FxHashMap::default(),
        }
    }

    fn len(&self) -> usize {
        self.entries.len()
    }

    fn mark_keys_dirty(&mut self) {
        self.dirty = true;
    }

    fn iter<'a>(
        &'a self,
    ) -> impl Iterator<Item = (&'a GeneralizedTableau<T, I, C>, &'a <T as Config>::Coeff)>
    where
        T: 'a,
        I: 'a,
        C: 'a,
    {
        self.entries.iter().map(|(t, c)| (t, c))
    }

    fn for_each_mut<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut <T as Config>::Coeff),
    {
        self.entries.iter_mut().for_each(|(tab, c)| f(tab, c));
    }

    fn for_each_mut_with_keys<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut <T as Config>::Coeff, u64, u64),
    {
        self.rebuild_fingerprints_if_dirty();
        for (i, (tab, c)) in self.entries.iter_mut().enumerate() {
            f(tab, c, self.word_fingerprints[i], self.phase_loss_hashes[i]);
        }
    }

    fn insert_or_merge_batch(
        &mut self,
        branches: Vec<(GeneralizedTableau<T, I, C>, <T as Config>::Coeff, u64, u64)>,
        cutoff: &<T as Config>::Coeff,
    ) -> bool {
        self.rebuild_fingerprints_if_dirty();

        // Build a fingerprint index over the current entries so each branch
        // lookup is O(1) expected instead of O(M) linear scan. Per-entry
        // fingerprints are cached on `self.fingerprints` and reused across
        // consecutive noise calls; Clifford/T gates clear the cache.
        let mut fp_index: FxHashMap<u64, Vec<usize>> =
            FxHashMap::with_capacity_and_hasher(self.entries.len(), Default::default());
        for i in 0..self.entries.len() {
            let fp = self.fingerprints[i];
            fp_index.entry(fp).or_default().push(i);
        }

        let mut needs_renormalize = false;
        for (tab, p, word_fp, phase_loss) in branches {
            let dropped_any =
                self.insert_or_merge_entry(tab, p, word_fp, phase_loss, &mut fp_index, cutoff);
            needs_renormalize |= dropped_any;
        }

        needs_renormalize
    }

    fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(&GeneralizedTableau<T, I, C>, &<T as Config>::Coeff) -> bool,
    {
        // Walk entries + fingerprints in lockstep so they stay aligned.
        // Order-preserving, O(N), no allocation.
        let mut write = 0;
        for read in 0..self.entries.len() {
            let (tab, c) = &self.entries[read];
            if f(tab, c) {
                if read != write {
                    self.entries.swap(read, write);
                    self.fingerprints.swap(read, write);
                    self.word_fingerprints.swap(read, write);
                    self.phase_loss_hashes.swap(read, write);
                }
                write += 1;
            }
        }
        self.entries.truncate(write);
        self.fingerprints.truncate(write);
        self.word_fingerprints.truncate(write);
        self.phase_loss_hashes.truncate(write);
    }

    fn reset_loss_and_merge(&mut self, addr0: usize) -> bool {
        self.rebuild_fingerprints_if_dirty();

        let delta = loss_mask(addr0);
        let mut changed_indices = Vec::new();
        let mut target_fingerprints: FxHashSet<u64> = FxHashSet::default();
        for (i, (tab, _)) in self.entries.iter().enumerate() {
            if tab.is_lost[addr0] {
                changed_indices.push(i);
                target_fingerprints.insert(self.fingerprints[i] ^ delta);
            }
        }
        if changed_indices.is_empty() {
            return false;
        }

        let mut changed = vec![false; self.entries.len()];
        for &i in &changed_indices {
            changed[i] = true;
        }

        let mut fp_index: FxHashMap<u64, SmallVec<[usize; 1]>> =
            FxHashMap::with_capacity_and_hasher(target_fingerprints.len(), Default::default());
        for i in 0..self.entries.len() {
            if changed[i] {
                continue;
            }

            let fp = self.fingerprints[i];
            if target_fingerprints.contains(&fp) {
                fp_index.entry(fp).or_default().push(i);
            }
        }

        let mut remove = vec![false; self.entries.len()];
        let mut merged_any = false;
        for i in changed_indices {
            self.entries[i].0.is_lost[addr0] = false;
            self.phase_loss_hashes[i] ^= delta;
            self.fingerprints[i] ^= delta;
            let fp = self.fingerprints[i];

            let found = {
                let tab = &self.entries[i].0;
                let scratch = &mut self.scratch;
                fp_index.get(&fp).and_then(|candidates| {
                    candidates.iter().copied().find(|&candidate| {
                        structurally_equal(&self.entries[candidate].0, tab, scratch)
                    })
                })
            };

            match found {
                Some(candidate) => {
                    self.entries[candidate].1 =
                        self.entries[candidate].1.clone() + self.entries[i].1.clone();
                    remove[i] = true;
                    merged_any = true;
                }
                None => fp_index.entry(fp).or_default().push(i),
            }
        }

        if merged_any {
            let mut write = 0;
            for read in 0..self.entries.len() {
                if !remove[read] {
                    if read != write {
                        self.entries.swap(read, write);
                        self.fingerprints.swap(read, write);
                        self.word_fingerprints.swap(read, write);
                        self.phase_loss_hashes.swap(read, write);
                    }
                    write += 1;
                }
            }
            self.entries.truncate(write);
            self.fingerprints.truncate(write);
            self.word_fingerprints.truncate(write);
            self.phase_loss_hashes.truncate(write);
        }

        merged_any
    }
}
