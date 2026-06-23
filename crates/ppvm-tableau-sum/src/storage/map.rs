// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::ops::AddAssign;

use fxhash::FxHashMap;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use ppvm_traits::config::Config;
use smallvec::SmallVec;

use crate::storage::{
    Branch, BranchMutation, EntryStore, apply_branch_mutation, fingerprint, phase_loss_hash,
    structurally_equal, word_fingerprint,
};
use bitvec::view::BitView;
use num::PrimInt;
use ppvm_traits::traits::Clifford;

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
    GeneralizedTableau<T, I, C>: Clifford,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
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

    fn insert_or_merge_mutated_branches(
        &mut self,
        branches: Vec<(usize, BranchMutation, <T as Config>::Coeff, u64, u64)>,
        cutoff: &<T as Config>::Coeff,
    ) -> bool {
        self.rebuild_if_dirty();

        // Materialize parents in the SAME order as for_each_mut_with_keys so
        // parent_idx aligns. Correctness-only path: no clone savings here.
        let parents: Vec<_> = self
            .buckets
            .values()
            .flat_map(|v| v.iter())
            .map(|(t, _)| t.clone())
            .collect();

        let real: Vec<Branch<T, I, C>> = branches
            .into_iter()
            .map(|(parent_idx, mutation, p, word_fp, phase_loss)| {
                let mut tab = parents[parent_idx].clone();
                apply_branch_mutation(&mut tab, mutation);
                (tab, p, word_fp, phase_loss)
            })
            .collect();

        self.insert_or_merge_batch(real, cutoff)
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

    fn drain_where<F>(
        &mut self,
        mut pred: F,
    ) -> Vec<(GeneralizedTableau<T, I, C>, <T as Config>::Coeff, u64, u64)>
    where
        F: FnMut(&GeneralizedTableau<T, I, C>) -> bool,
    {
        self.rebuild_if_dirty();

        let mut drained = Vec::new();
        // Buckets store entries under the full fingerprint, so the split is
        // (fp, 0). XORing a delta into either component during the caller's
        // mutation produces the right new bucket key in insert_or_merge_batch.
        self.buckets.retain(|&fp, bucket| {
            let mut i = 0;
            while i < bucket.len() {
                if pred(&bucket[i].0) {
                    let (tab, p) = bucket.swap_remove(i);
                    drained.push((tab, p, fp, 0));
                } else {
                    i += 1;
                }
            }
            !bucket.is_empty()
        });
        drained
    }
}
