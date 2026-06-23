// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashSet;

use ppvm_traits::config::Config;
use ppvm_traits::traits::*;

/// A sparse formal sum `Σ cᵢ Pᵢ` of Pauli strings.
///
/// The central data type of the Pauli-propagation backend. Keys are
/// [`PauliWord`](ppvm_pauli_word::word::PauliWord)-shaped Pauli strings; values are
/// numeric coefficients. Generic over a [`Config`] that fixes the
/// concrete storage / hasher / coefficient / strategy.
///
/// Internally `PauliSum` holds *two* maps — a primary and an auxiliary —
/// and swaps between them during gate propagation to avoid repeated
/// allocations. The [`PauliSum::data`] / [`PauliSum::aux`] accessors
/// expose the current orientation; most callers will only ever touch
/// the primary side via the high-level gate / measurement traits.
///
/// # Examples
///
/// Heisenberg-picture propagation of `ZZ` through the GHZ circuit:
///
/// ```
/// use ppvm_pauli_sum::prelude::*;
///
/// let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
///     PauliSum::builder().n_qubits(2).build();
/// state += ("ZZ", 1.0);
///
/// // Circuit: H(0); CNOT(0, 1) — apply in reverse for Heisenberg propagation.
/// state.cnot(0, 1);
/// state.h(0);
///
/// // ZZ → IZ under the GHZ circuit, with coefficient 1.0.
/// assert_eq!(state.len(), 1);
/// ```
#[derive(Clone)]
pub struct PauliSum<T: Config> {
    map: (T::Map, T::Map),
    aux: bool,
    n_qubits: usize,
    capacity: usize,
    strategy: T::Strategy,
    /// Keep-set: strings [`PauliSum::truncate`] must always re-insert
    /// after the strategy runs. Empty by default.
    preserve_strings: HashSet<T::PauliWordType>,
    /// Reusable scratch buffer for [`map_insert`](Self::map_insert): produced
    /// entries are pushed here (no hashing) then merged into the primary map,
    /// avoiding a second hashmap probe per entry.
    scratch: Vec<(T::PauliWordType, T::Coeff)>,
}

#[bon::bon]
impl<T: Config> PauliSum<T> {
    /// create a new empty PauliSum with given number of qubits.
    ///
    /// One can optionally set
    /// - the strategy for truncation, initialization etc.
    /// - the capacity of the internal maps, default is strategy.capacity(n_qubits)
    /// - a set of `preserve_strings` that [`truncate`](Self::truncate) must
    ///   never drop, on top of whatever the strategy decides.
    #[builder]
    pub fn new(
        /// number of qubits
        n_qubits: usize,
        /// strategy for truncation, initialization etc.
        #[builder(default = T::Strategy::default())]
        strategy: T::Strategy,
        /// capacity of the internal maps, default is strategy.capacity(n_qubits)
        #[builder(default = strategy.capacity(n_qubits))]
        capacity: usize,
        /// Pauli strings that truncate must always keep. Empty by default.
        #[builder(default)]
        preserve_strings: HashSet<T::PauliWordType>,
    ) -> Self {
        Self {
            map: (
                T::Map::with_capacity(capacity),
                T::Map::with_capacity(capacity),
            ),
            aux: false,
            n_qubits,
            capacity,
            strategy,
            preserve_strings,
            scratch: Vec::new(),
        }
    }
}

impl<T: Config> PauliSum<T> {
    /// Number of qubits the sum is defined over.
    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Capacity (in entries) reserved in the underlying maps.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Reference to the primary map — the one that currently holds the
    /// "active" entries.
    #[inline(always)]
    pub fn data(&self) -> &T::Map {
        if self.aux { &self.map.1 } else { &self.map.0 }
    }

    /// Mutable reference to the primary map.
    #[inline(always)]
    pub fn data_mut(&mut self) -> &mut T::Map {
        if self.aux {
            &mut self.map.1
        } else {
            &mut self.map.0
        }
    }

    /// Reference to the auxiliary map — the scratch buffer gates write into.
    #[inline(always)]
    pub fn aux(&self) -> &T::Map {
        if self.aux { &self.map.0 } else { &self.map.1 }
    }

    /// Mutable reference to the auxiliary map.
    #[inline(always)]
    pub fn aux_mut(&mut self) -> &mut T::Map {
        if self.aux {
            &mut self.map.0
        } else {
            &mut self.map.1
        }
    }

    /// Both maps as a `(primary, auxiliary)` tuple.
    #[inline(always)]
    pub fn data_aux(&self) -> (&T::Map, &T::Map) {
        if self.aux {
            (&self.map.1, &self.map.0)
        } else {
            (&self.map.0, &self.map.1)
        }
    }

    /// Both maps mutably as a `(primary, auxiliary)` tuple.
    #[inline(always)]
    pub fn data_aux_mut(&mut self) -> (&mut T::Map, &mut T::Map) {
        if self.aux {
            (&mut self.map.1, &mut self.map.0)
        } else {
            (&mut self.map.0, &mut self.map.1)
        }
    }

    /// Swap the primary and auxiliary maps without touching their contents.
    #[inline(always)]
    pub fn swap(&mut self) {
        self.aux = !self.aux;
    }

    /// Number of entries in the primary map.
    pub fn len(&self) -> usize {
        self.data().len()
    }

    /// `true` if the primary map has no entries.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `true` if `(key, value)` is present in the primary map.
    pub fn contains(&self, key: &<T as Config>::PauliWordType, value: &T::Coeff) -> bool {
        self.data().contains(key, value)
    }

    /// combine entries with the same key
    /// in either data or aux. The combined entries are stored in `.data()`.
    pub fn consume(&mut self) {
        let (data, aux) = self.data_aux_mut();
        if aux.len() > data.len() {
            aux.consume(data);
            self.swap();
        } else {
            data.consume(aux);
        }
    }

    /// scale all coefficients by a function of the corresponding PauliWord key.
    pub fn scale<F>(&mut self, f: F)
    where
        F: Fn(&<T as Config>::PauliWordType, &mut T::Coeff) + Sync + Send,
    {
        self.data_mut().scale(f);
    }

    /// modify in place existing entries and insert some new entries
    /// if `f` return Some((k,v)) for an existing entry (k0,v0), then
    /// the existing entry is modified by `f` and a new entry (k,v) is added.
    /// if `f` return None, then the existing entry is only modified.
    /// finally, all entries are combined assuming unique keys.
    pub fn map_insert<F>(&mut self, f: F)
    where
        F: Fn(
                &<T as Config>::PauliWordType,
                &mut T::Coeff,
            ) -> Option<(<T as Config>::PauliWordType, T::Coeff)>
            + Sync
            + Send,
    {
        // Collect produced terms into a reusable scratch Vec (cheap pushes,
        // no hashing), then merge them into the primary map in one pass. This
        // avoids the extra hashmap probe + insert the old aux-map buffer paid
        // for every produced term, and makes the per-call clear O(1).
        let mut scratch = std::mem::take(&mut self.scratch);
        scratch.clear();
        self.data_mut().map_insert_vec(&mut scratch, f);

        let data = self.data_mut();
        for (k, v) in scratch.drain(..) {
            data.add_assign(k, v);
        }
        self.scratch = scratch;
    }

    /// modify in place existing entries and insert some new entries
    /// if `f` returns Some(Vec<(k,v)>) for an existing entry (k0,v0), then
    /// the existing entry is modified by `f` and new entries contained in the
    /// vector are added.
    /// If `f` returns None, then the existing entry is only modified.
    /// finally, all entries are combined assuming unique keys.
    pub fn map_insert_multiple<F>(&mut self, f: F)
    where
        F: Fn(
                &<T as Config>::PauliWordType,
                &mut T::Coeff,
            ) -> Option<Vec<(<T as Config>::PauliWordType, T::Coeff)>>
            + Sync
            + Send,
    {
        let (data, aux) = self.data_aux_mut();
        aux.clear();
        data.map_insert_multiple(aux, f);
        self.consume();
    }

    /// apply a function to each entry (k,v) and store the results in aux.
    /// finally, swap data and aux.
    ///
    /// This assumes the function `f` returns a different Pauli string `k'`
    /// from the input `k`, otherwise use `scale` to modify coefficients in place
    /// of the same entry.
    pub fn map_add<F>(&mut self, f: F)
    where
        F: Fn(&<T as Config>::PauliWordType, &T::Coeff) -> (<T as Config>::PauliWordType, T::Coeff)
            + Sync
            + Send,
    {
        let (data, aux) = self.data_aux_mut();
        aux.clear();
        data.map_add_assign(aux, f);
        self.swap();
    }

    /// Apply the configured truncation [`Strategy`](ppvm_traits::traits::Strategy)
    /// to the primary map, dropping entries that fall outside its policy.
    ///
    /// If `preserve_strings` is non-empty, any of those Pauli strings
    /// that the strategy would have dropped are re-inserted afterwards
    /// with their pre-truncate coefficient. The mechanism composes with
    /// any [`Strategy`] (coefficient-magnitude, max-weight, combinations
    /// — anything) because the strategy runs unchanged in the middle.
    pub fn truncate(&mut self) {
        // Hot path: empty preserve set → just run the strategy.
        if self.preserve_strings.is_empty() {
            let strategy = self.strategy;
            strategy.truncate(self.data_mut());
            return;
        }

        // Snapshot the current coefficients of preserved keys. We piggy-
        // back on `retain` (which always returns true here, so it's a
        // pure scan) to walk `(k, v)` pairs without needing a separate
        // `get`/`iter` route through the `ACMap` traits.
        let preserve = self.preserve_strings.clone();
        let mut saved: Vec<(T::PauliWordType, T::Coeff)> = Vec::new();
        self.data_mut().retain(|k, v| {
            if preserve.contains(k) {
                saved.push((k.clone(), v.clone()));
            }
            true
        });

        // Run the configured strategy verbatim.
        let strategy = self.strategy;
        strategy.truncate(self.data_mut());

        // Restore any preserved entry the strategy dropped. We use
        // `add_assign` because the trait does not expose a plain
        // "insert if absent"; on a missing key it inserts, which is
        // what we want here (the guarded `contains_with` keeps us from
        // accidentally summing into a kept entry).
        let data = self.data_mut();
        for (k, v) in saved {
            if !data.contains_with(&k, |_| true) {
                data.add_assign(k, v);
            }
        }
    }

    /// Read-only access to the active preserve set.
    pub fn preserve_strings(&self) -> &HashSet<T::PauliWordType> {
        &self.preserve_strings
    }
}

impl<'a, T: Config> PauliSum<T>
where
    T::Map: ACMapIter<'a>,
{
    /// Iterate over the primary map's entries.
    pub fn iter(&'a self) -> <<T as Config>::Map as ACMapIter<'a>>::Iter {
        self.data().iter()
    }
}

impl<T: Config> IntoIterator for PauliSum<T>
where
    T::Map: IntoIterator,
{
    type Item = <T::Map as IntoIterator>::Item;
    type IntoIter = <T::Map as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        if self.aux {
            self.map.1.into_iter()
        } else {
            self.map.0.into_iter()
        }
    }
}

impl<T: Config> PartialEq for PauliSum<T>
where
    T::Map: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.n_qubits == other.n_qubits && self.data() == other.data()
    }
}

impl<T: Config> Extend<(<T as Config>::PauliWordType, T::Coeff)> for PauliSum<T>
where
    T::Map: Extend<(<T as Config>::PauliWordType, T::Coeff)>,
{
    fn extend<I: IntoIterator<Item = (<T as Config>::PauliWordType, T::Coeff)>>(
        &mut self,
        iter: I,
    ) {
        self.data_mut().extend(iter);
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_yaml_snapshot;

    use super::*;
    use crate::config::fxhash::ByteF64;
    use ppvm_pauli_word::word::PauliWord;

    #[test]
    fn test_pauli_sum_creation() {
        let word = PauliWord::<[u8; 2]>::new(4);
        let mut sum: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(word.n_qubits()).build();
        assert!(sum.data().is_empty());
        sum += "IIII";
        assert!(!sum.data().is_empty());
        assert_yaml_snapshot!(sum.to_string());
        sum += ("IIII", 2.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += ("XIII", 1.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += ("XIII", 2.0);
        assert_yaml_snapshot!(sum.to_string());
        sum += "IYII";
        assert_yaml_snapshot!(sum.to_string());
        assert!(sum.contains(&PauliWord::from("IIII"), &3.0));
        assert!(sum.contains(&PauliWord::from("XIII"), &3.0));
        assert!(sum.contains(&PauliWord::from("IYII"), &1.0));
    }

    #[test]
    fn test_pauli_sum_top_bottom() {
        let mut sum: PauliSum<ByteF64<2>> = PauliSum::builder().n_qubits(4).build();
        assert!(sum.is_empty());
        sum += ("IIII", 1.0);
        assert!(!sum.is_empty());
        sum += ("IIII", 1.0);
    }
}
