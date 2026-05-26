use std::ops::AddAssign;

use fxhash::FxHashMap;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};

use crate::storage::{EntryStore, fingerprint, structurally_equal};

#[derive(Clone)]
pub struct VecStorage<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> {
    pub entries: Vec<(GeneralizedTableau<T, I, C>, T::Coeff)>,
    pub fingerprints: Vec<u64>,
    pub dirty: bool,
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
        fp_index: &mut FxHashMap<u64, Vec<usize>>,
        sum_cutoff: &T::Coeff,
    ) -> bool {
        // New branches always come from a freshly-mutated fork (X/Y/Z or
        // is_lost flip applied just before this call), so their cached
        // fingerprint would always be stale — compute fresh.
        let fp = fingerprint(&tab);

        // Only run the full equality check on entries whose fingerprint matches.
        let mut found: Option<usize> = None;
        if let Some(candidates) = fp_index.get(&fp) {
            for &i in candidates {
                if structurally_equal(&self.entries[i].0, &tab) {
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
                    fp_index.entry(fp).or_default().push(new_idx);
                } else {
                    needs_normalize = true;
                }
            }
        }

        needs_normalize
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
            dirty: false,
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

    fn insert_or_merge_batch(
        &mut self,
        branches: Vec<(GeneralizedTableau<T, I, C>, <T as Config>::Coeff)>,
        cutoff: &<T as Config>::Coeff,
    ) -> bool {
        if self.dirty {
            self.fingerprints.clear();
            self.fingerprints
                .extend(self.entries.iter().map(|(t, _)| fingerprint(t)));
            self.dirty = false;
        }

        // Build a fingerprint index over the current entries so each branch
        // lookup is O(1) expected instead of O(M) linear scan. Per-entry
        // fingerprints are cached on `self.entry_fingerprints` and reused
        // across consecutive noise calls; Clifford/T gates clear the cache.
        let mut fp_index: FxHashMap<u64, Vec<usize>> =
            FxHashMap::with_capacity_and_hasher(self.entries.len(), Default::default());
        for i in 0..self.entries.len() {
            let fp = self.fingerprints[i];
            fp_index.entry(fp).or_default().push(i);
        }

        let mut needs_renormalize = false;
        for (tab, p) in branches {
            let dropped_any = self.insert_or_merge_entry(tab, p, &mut fp_index, cutoff);
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
                }
                write += 1;
            }
        }
        self.entries.truncate(write);
        self.fingerprints.truncate(write);
    }
}
