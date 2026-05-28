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
use smallvec::SmallVec;

use crate::storage::{
    EntryStore, fingerprint, phase_loss_hash, structurally_equal, word_fingerprint,
};

type Bucket<T, I, C> = SmallVec<[(GeneralizedTableau<T, I, C>, <T as Config>::Coeff); 1]>;

#[derive(Clone)]
pub struct MapStorage<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> {
    pub buckets: FxHashMap<u64, Bucket<T, I, C>>,
    pub dirty: bool,
    /// Reusable scratch buffer for `structurally_equal`'s coefficient lookup
    /// map. Cleared and refilled per call; keeps its capacity across calls.
    scratch: FxHashMap<I, Complex<T::Coeff>>,
}

impl<T, I, C> MapStorage<T, I, C>
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
    /// Re-bucket every entry under its fresh fingerprint. Clifford/T gates
    /// mutate tableau data in place, so after `mark_keys_dirty` the existing
    /// keys are wrong. Iterates the old map directly to avoid materializing
    /// the entries into a temporary Vec.
    fn rebuild_if_dirty(&mut self) {
        if !self.dirty {
            return;
        }
        let entry_count = self.buckets.values().map(|v| v.len()).sum();
        let old = std::mem::replace(
            &mut self.buckets,
            FxHashMap::with_capacity_and_hasher(entry_count, Default::default()),
        );
        for (_, bucket) in old {
            for (tab, p) in bucket {
                let fp = fingerprint(&tab);
                self.buckets.entry(fp).or_default().push((tab, p));
            }
        }
        self.dirty = false;
    }
}

impl<T, I, C> EntryStore<T, I, C> for MapStorage<T, I, C>
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
            buckets: FxHashMap::with_capacity_and_hasher(cap, Default::default()),
            dirty: false,
            scratch: FxHashMap::default(),
        }
    }

    fn len(&self) -> usize {
        self.buckets.values().map(|v| v.len()).sum()
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
        self.buckets
            .values()
            .flat_map(|v| v.iter())
            .map(|(t, c)| (t, c))
    }

    fn for_each_mut<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut <T as Config>::Coeff),
    {
        for v in self.buckets.values_mut() {
            for (tab, c) in v.iter_mut() {
                f(tab, c);
            }
        }
    }

    fn for_each_mut_with_keys<F>(&mut self, mut f: F)
    where
        F: FnMut(&mut GeneralizedTableau<T, I, C>, &mut <T as Config>::Coeff, u64, u64),
    {
        self.rebuild_if_dirty();
        for v in self.buckets.values_mut() {
            for (tab, c) in v.iter_mut() {
                let word_fp = word_fingerprint(tab);
                let phase_loss = phase_loss_hash(tab);
                f(tab, c, word_fp, phase_loss);
            }
        }
    }

    fn insert_or_merge_batch(
        &mut self,
        branches: Vec<(GeneralizedTableau<T, I, C>, <T as Config>::Coeff, u64, u64)>,
        cutoff: &<T as Config>::Coeff,
    ) -> bool {
        self.rebuild_if_dirty();

        let mut needs_renormalize = false;
        for (tab, p, word_fp, phase_loss) in branches {
            // The branch carries both fingerprint components, so the full
            // fingerprint is their XOR — no walk over the tableau.
            let fp = word_fp ^ phase_loss;
            let bucket = self.buckets.entry(fp).or_default();

            let mut found: Option<usize> = None;
            for (i, (existing, _)) in bucket.iter().enumerate() {
                if structurally_equal(existing, &tab, &mut self.scratch) {
                    found = Some(i);
                    break;
                }
            }

            match found {
                Some(i) => {
                    let p0 = &bucket[i].1;
                    bucket[i].1 = p0.clone() + p;
                }
                None => {
                    if &p > cutoff {
                        bucket.push((tab, p));
                    } else {
                        // Drop the branch. If we just created the bucket via
                        // `entry().or_default()`, remove it so iter() stays
                        // clean and len() isn't skewed by empty slots.
                        if bucket.is_empty() {
                            self.buckets.remove(&fp);
                        }
                        needs_renormalize = true;
                    }
                }
            }
        }

        needs_renormalize
    }

    fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(&GeneralizedTableau<T, I, C>, &<T as Config>::Coeff) -> bool,
    {
        self.buckets.retain(|_, v| {
            v.retain(|(tab, c)| f(tab, c));
            !v.is_empty()
        });
    }

    fn merge_equal_entries(&mut self) -> bool {
        let entry_count = self.len();
        if entry_count < 2 {
            self.rebuild_if_dirty();
            return false;
        }

        let old = std::mem::replace(
            &mut self.buckets,
            FxHashMap::with_capacity_and_hasher(entry_count, Default::default()),
        );
        let mut merged_any = false;

        for (_, bucket) in old {
            for (tab, p) in bucket {
                let fp = fingerprint(&tab);
                let bucket = self.buckets.entry(fp).or_default();

                let mut found: Option<usize> = None;
                for (i, (existing, _)) in bucket.iter().enumerate() {
                    if structurally_equal(existing, &tab, &mut self.scratch) {
                        found = Some(i);
                        break;
                    }
                }

                match found {
                    Some(i) => {
                        bucket[i].1 = bucket[i].1.clone() + p;
                        merged_any = true;
                    }
                    None => bucket.push((tab, p)),
                }
            }
        }

        self.dirty = false;
        merged_any
    }
}
