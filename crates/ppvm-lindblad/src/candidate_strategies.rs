// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Bounded-memory alternatives to chunk-capped signed leakage accumulation.

use crate::algebra::{anti_commutes, comm_product};
use crate::direct_expm::ActionScratch;
use crate::indexed_candidates::IndexedCandidates;
use crate::word::word_hash;
use crate::{CandidateBuffer, CandidateProposer, Error, LindbladSpec, StreamingProposer, Word};
use fxhash::FxHashSet;
use hashbrown::HashTable;
use rayon::prelude::*;

fn source_order(order: &mut Vec<usize>, coefficients: &[f64]) {
    order.clear();
    order.extend(0..coefficients.len());
    order.sort_unstable_by(|&a, &b| {
        coefficients[b]
            .abs()
            .total_cmp(&coefficients[a].abs())
            .then_with(|| a.cmp(&b))
    });
}

#[derive(Default)]
struct Positions {
    table: HashTable<u32>,
}

fn positioned_word<'a, const C: usize>(
    buffer: &'a CandidateBuffer<'_, C>,
    position: u32,
) -> &'a Word<C> {
    let position = position as usize;
    if position < buffer.basis().len() {
        &buffer.basis()[position]
    } else {
        &buffer.candidates()[position - buffer.basis().len()]
    }
}

impl Positions {
    fn reset<const C: usize>(&mut self, buffer: &CandidateBuffer<'_, C>) {
        self.table.clear();
        for (index, word) in buffer.basis().iter().enumerate() {
            self.table.insert_unique(
                word_hash(word),
                u32::try_from(index).expect("basis exceeds u32 positions"),
                |&index| word_hash(positioned_word(buffer, index)),
            );
        }
    }

    fn get<const C: usize>(&self, buffer: &CandidateBuffer<'_, C>, word: &Word<C>) -> Option<u32> {
        self.table
            .find(word_hash(word), |&index| {
                positioned_word(buffer, index) == word
            })
            .copied()
    }

    fn insert_candidate<const C: usize>(
        &mut self,
        buffer: &CandidateBuffer<'_, C>,
        candidate: usize,
    ) {
        let position = u32::try_from(buffer.basis().len() + candidate)
            .expect("candidate positions exceed u32");
        self.table.insert_unique(
            word_hash(&buffer.candidates()[candidate]),
            position,
            |&index| word_hash(positioned_word(buffer, index)),
        );
    }
}

#[derive(Default)]
struct ScoreHeap {
    limit: usize,
    capacity: usize,
    len: usize,
}

impl ScoreHeap {
    fn ensure(&mut self, numeric: &mut Vec<f64>, required: usize) {
        if required <= self.capacity {
            return;
        }
        let previous = self.capacity;
        self.capacity = required
            .max(previous.saturating_mul(2))
            .max(8)
            .min(self.limit);
        numeric.resize(2 * self.capacity, 0.0);
        if previous != 0 {
            numeric.copy_within(previous..2 * previous, self.capacity);
        }
    }
    // The packed region stores heap slots in low 32 bits and inverse positions
    // in high 32 bits. Floating-point arithmetic only reads the score prefix.
    fn slot(&self, numeric: &[f64], rank: usize) -> usize {
        numeric[self.capacity + rank].to_bits() as u32 as usize
    }
    fn rank(&self, numeric: &[f64], slot: usize) -> usize {
        (numeric[self.capacity + slot].to_bits() >> 32) as usize
    }
    fn set_slot(&self, numeric: &mut [f64], rank: usize, slot: usize) {
        let old = numeric[self.capacity + rank].to_bits();
        numeric[self.capacity + rank] = f64::from_bits(
            (old & 0xffffffff00000000) | u32::try_from(slot).expect("heap slot exceeds u32") as u64,
        );
    }
    fn set_rank(&self, numeric: &mut [f64], slot: usize, rank: usize) {
        let old = numeric[self.capacity + slot].to_bits();
        numeric[self.capacity + slot] = f64::from_bits(
            (old & 0xffffffff)
                | ((u32::try_from(rank).expect("heap rank exceeds u32") as u64) << 32),
        );
    }
    fn worse(numeric: &[f64], a: usize, b: usize) -> bool {
        numeric[a].abs() < numeric[b].abs() || (numeric[a].abs() == numeric[b].abs() && a > b)
    }
    fn swap(&self, numeric: &mut [f64], a: usize, b: usize) {
        let sa = self.slot(numeric, a);
        let sb = self.slot(numeric, b);
        self.set_slot(numeric, a, sb);
        self.set_slot(numeric, b, sa);
        self.set_rank(numeric, sa, b);
        self.set_rank(numeric, sb, a);
    }
    fn insert(&mut self, numeric: &mut Vec<f64>, slot: usize, value: f64) {
        self.ensure(numeric, slot + 1);
        numeric[slot] = value;
        let mut rank = self.len;
        self.len += 1;
        self.set_slot(numeric, rank, slot);
        self.set_rank(numeric, slot, rank);
        while rank > 0 {
            let parent = (rank - 1) / 2;
            if !Self::worse(
                numeric,
                self.slot(numeric, rank),
                self.slot(numeric, parent),
            ) {
                break;
            }
            self.swap(numeric, rank, parent);
            rank = parent;
        }
    }
    fn update(&self, numeric: &mut [f64], slot: usize, value: f64) {
        numeric[slot] = value;
        let mut rank = self.rank(numeric, slot);
        loop {
            let left = 2 * rank + 1;
            if left >= self.len {
                break;
            }
            let right = left + 1;
            let child = if right < self.len
                && Self::worse(numeric, self.slot(numeric, right), self.slot(numeric, left))
            {
                right
            } else {
                left
            };
            if !Self::worse(numeric, self.slot(numeric, child), self.slot(numeric, rank)) {
                break;
            }
            self.swap(numeric, rank, child);
            rank = child;
        }
    }
    fn minimum(&self, numeric: &[f64]) -> (usize, f64) {
        let slot = self.slot(numeric, 0);
        (slot, numeric[slot].abs())
    }
}

/// Select distinct external strings by their largest individual `|h_i c_j|`.
///
/// A bounded minimum heap keeps the largest individual contribution per destination.
/// Sorted-product upper bounds skip pairs below its cutoff. Duplicate amplitudes
/// are never summed. Equal scores have deterministic, unspecified tie order.
/// Pauli-noise jumps are diagonal and create no external candidates. General
/// dissipators use the streaming proposer to retain support for those models.
#[derive(Default)]
pub struct GreedyProposer<const C: usize = { crate::W_CHUNKS }> {
    h: Vec<(Word<C>, f64)>,
    order: Vec<usize>,
    heap: ScoreHeap,
    positions: Positions,
    protected: FxHashSet<Word<C>>,
    fallback: StreamingProposer<C>,
    use_fallback: bool,
}

impl<const C: usize> GreedyProposer<C> {
    fn prepare(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &CandidateBuffer<'_, C>,
        protected: &[Word<C>],
    ) {
        self.heap.capacity = 0;
        self.heap.len = 0;
        self.h.clear();
        self.h
            .extend(model.hamiltonian_terms().filter(|(_, c)| *c != 0.0));
        self.h.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
        source_order(&mut self.order, buffer.coefficients());
        self.positions.reset(buffer);
        self.protected.clear();
        self.protected.extend(protected.iter().copied());
    }

    fn take_range(
        &mut self,
        buffer: &mut CandidateBuffer<'_, C>,
        start: usize,
        end: usize,
        budget: usize,
    ) {
        if budget == 0 || start == end {
            return;
        }
        self.heap.len = 0;
        let base = buffer.candidates().len();
        self.heap.limit = base.saturating_add(budget);
        for source_rank in start..end {
            let source = self.order[source_rank];
            let coefficient = buffer.coefficients()[source];
            for &(h_word, h_coefficient) in &self.h {
                let bound = 2.0 * h_coefficient.abs() * coefficient.abs();
                if self.heap.len == budget && bound <= self.heap.minimum(buffer.numeric).1 {
                    break;
                }
                if bound == 0.0 {
                    break;
                }
                if !anti_commutes(&h_word, &buffer.basis()[source]) {
                    continue;
                }
                let (word, sign) = comm_product(&h_word, &buffer.basis()[source]);
                if self.protected.contains(&word) {
                    continue;
                }
                let rate = h_coefficient * sign * coefficient;
                if let Some(position) = self.positions.get(buffer, &word) {
                    let position = position as usize;
                    if position < buffer.basis().len() + base {
                        continue;
                    }
                    let slot = position - buffer.basis().len();
                    if rate.abs() > buffer.numeric[slot].abs() {
                        self.heap.update(buffer.numeric, slot, rate);
                    }
                } else if self.heap.len < budget {
                    let slot = buffer.candidates().len();
                    buffer.push(word);
                    self.heap.insert(buffer.numeric, slot, rate);
                    self.positions.insert_candidate(buffer, slot);
                } else {
                    let slot = self.heap.minimum(buffer.numeric).0;
                    let position = u32::try_from(buffer.basis().len() + slot)
                        .expect("candidate position exceeds u32");
                    self.positions
                        .table
                        .find_entry(word_hash(&buffer.candidates()[slot]), |&p| p == position)
                        .expect("heap destination absent from position table")
                        .remove();
                    buffer.candidates_mut()[slot] = word;
                    self.positions.insert_candidate(buffer, slot);
                    self.heap.update(buffer.numeric, slot, rate);
                }
            }
            if let Some((_, h_coefficient)) = self.h.first()
                && self.heap.len == budget
                && source_rank + 1 < end
                && 2.0
                    * h_coefficient.abs()
                    * buffer.coefficients()[self.order[source_rank + 1]].abs()
                    <= self.heap.minimum(buffer.numeric).1
            {
                break;
            }
        }
    }
}

fn compact_threshold<const C: usize>(buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
    let mut write = 0;
    for read in 0..buffer.candidates().len() {
        if buffer.numeric[read].abs() > threshold {
            buffer
                .words
                .swap(buffer.prefix + write, buffer.prefix + read);
            buffer.numeric[write] = buffer.numeric[read];
            write += 1;
        }
    }
    write
}

impl<const C: usize> CandidateProposer<C> for GreedyProposer<C> {
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error> {
        self.use_fallback = !model.supports_row_gather();
        if self.use_fallback {
            return self.fallback.propose(model, buffer, protected, maximum);
        }
        let room = maximum.saturating_sub(buffer.basis().len());
        if room == 0 {
            return Ok(());
        }
        self.prepare(model, buffer, protected);
        self.take_range(buffer, 0, self.order.len(), room);
        Ok(())
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
        if self.use_fallback {
            self.fallback.select(buffer, threshold)
        } else {
            buffer.numeric.truncate(buffer.candidates().len());
            compact_threshold(buffer, threshold)
        }
    }
}

/// Divide the remaining candidate budget across remaining source batches.
///
/// Sources are sorted by coefficient magnitude. Each batch selects the largest
/// individual Hamiltonian/source contributions for distinct destinations. Words
/// already admitted by earlier batches are skipped without merging amplitudes.
pub struct BatchedProposer<const C: usize = { crate::W_CHUNKS }> {
    batches: usize,
    greedy: GreedyProposer<C>,
}

impl<const C: usize> BatchedProposer<C> {
    /// Use `batches` source partitions, clamped to at least one.
    pub fn new(batches: usize) -> Self {
        Self {
            batches: batches.max(1),
            greedy: GreedyProposer::default(),
        }
    }
}
impl<const C: usize> Default for BatchedProposer<C> {
    fn default() -> Self {
        Self::new(32)
    }
}
impl<const C: usize> CandidateProposer<C> for BatchedProposer<C> {
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error> {
        self.greedy.use_fallback = !model.supports_row_gather();
        if self.greedy.use_fallback {
            return self
                .greedy
                .fallback
                .propose(model, buffer, protected, maximum);
        }
        let room = maximum.saturating_sub(buffer.basis().len());
        if room == 0 || buffer.basis().is_empty() {
            return Ok(());
        }
        self.greedy.prepare(model, buffer, protected);
        let sources = self.greedy.order.len();
        let batches = self.batches.min(sources);
        for batch in 0..batches {
            let remaining = room - buffer.candidates().len();
            let budget = remaining / (batches - batch);
            self.greedy.take_range(
                buffer,
                sources * batch / batches,
                sources * (batch + 1) / batches,
                budget,
            );
        }
        Ok(())
    }
    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
        self.greedy.select(buffer, threshold)
    }
}

/// Build a greedy shortlist, compute its exact signed leakage, then rank it.
///
/// Exact scores use reverse Hamiltonian actions and a basis-position lookup,
/// including all contributions and cancellations for each shortlisted word.
/// Selection is exact within the shortlist, not over all external words.
pub struct RefinedGreedyProposer<const C: usize = { crate::W_CHUNKS }> {
    extra_percent: usize,
    room: usize,
    greedy: GreedyProposer<C>,
}
impl<const C: usize> RefinedGreedyProposer<C> {
    /// Oversample the admission budget by `extra_percent` before exact scoring.
    pub fn new(extra_percent: usize) -> Self {
        Self {
            extra_percent,
            room: 0,
            greedy: GreedyProposer::default(),
        }
    }
}
impl<const C: usize> Default for RefinedGreedyProposer<C> {
    fn default() -> Self {
        Self::new(50)
    }
}

fn select_largest<const C: usize>(
    buffer: &mut CandidateBuffer<'_, C>,
    limit: usize,
    threshold: f64,
) -> usize {
    let n = buffer.candidates().len();
    if n <= limit {
        return compact_threshold(buffer, threshold);
    }
    if limit == 0 {
        return 0;
    }
    buffer.numeric.resize(2 * n, 0.0);
    let (rates, keys) = buffer.numeric.split_at_mut(n);
    for (key, rate) in keys.iter_mut().zip(rates.iter()) {
        *key = rate.abs();
    }
    let cutoff = *keys.select_nth_unstable_by(n - limit, f64::total_cmp).1;
    let mut ties = limit - keys.iter().filter(|&&key| key > cutoff).count();
    buffer.numeric.truncate(n);
    let mut write = 0;
    for read in 0..n {
        let magnitude = buffer.numeric[read].abs();
        let keep = magnitude > cutoff || (magnitude == cutoff && ties > 0);
        if magnitude == cutoff && ties > 0 {
            ties -= 1;
        }
        if keep && magnitude > threshold {
            buffer
                .words
                .swap(buffer.prefix + write, buffer.prefix + read);
            buffer.numeric[write] = buffer.numeric[read];
            write += 1;
        }
    }
    write
}

impl<const C: usize> CandidateProposer<C> for RefinedGreedyProposer<C> {
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error> {
        self.room = maximum.saturating_sub(buffer.basis().len());
        self.greedy.use_fallback = !model.supports_row_gather();
        if self.greedy.use_fallback {
            return self
                .greedy
                .fallback
                .propose(model, buffer, protected, maximum);
        }
        if self.room == 0 {
            return Ok(());
        }
        self.greedy.prepare(model, buffer, protected);
        let shortlist = self
            .room
            .saturating_add(self.room.saturating_mul(self.extra_percent) / 100);
        self.greedy
            .take_range(buffer, 0, self.greedy.order.len(), shortlist);
        buffer.numeric.truncate(buffer.candidates().len());
        let (basis, words, coefficients, numeric) = buffer.scoring_parts();
        let n = basis.len();
        let positions = &self.greedy.positions;
        // The position table contains basis and shortlist entries. Only the
        // original basis has nonzero coefficients in the leakage calculation.
        let table = &positions.table;
        numeric.par_iter_mut().zip(words.par_iter()).for_each_init(
            ActionScratch::new,
            |scratch, (rate, word)| {
                scratch.fill(model, word, true);
                scratch.row.clear();
                for &(parent, value) in &scratch.terms {
                    if let Some(&index) = table.find(word_hash(&parent), |&index| {
                        (index as usize) < n && basis[index as usize] == parent
                    }) {
                        scratch.row.push((index, -value));
                    }
                }
                scratch.row.sort_unstable_by_key(|entry| entry.0);
                *rate = scratch
                    .row
                    .iter()
                    .map(|&(index, value)| coefficients[index as usize] * value)
                    .sum();
            },
        );
        Ok(())
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
        if self.greedy.use_fallback {
            self.greedy.fallback.select(buffer, threshold)
        } else {
            select_largest(buffer, self.room, threshold)
        }
    }
}

/// Merge signed leakage and prune only when the candidate pool exceeds headroom.
///
/// Small parallel source batches feed signed sums in source order. The default 25% headroom reduces whole-pool scans while fitting
/// N + 1.25 A candidate words into the usual 4B reservation when N=B and A=2B.
pub struct DeferredPruningProposer<const C: usize = { crate::W_CHUNKS }> {
    extra_percent: usize,
    indexed: IndexedCandidates,
    positions: Positions,
    protected: FxHashSet<Word<C>>,
    order: Vec<usize>,
    scratch: ActionScratch<C>,
    generation_batch: usize,
    outputs: Vec<Vec<(Word<C>, f64)>>,
}
impl<const C: usize> DeferredPruningProposer<C> {
    /// Allow `extra_percent` additional live candidates between pruning passes.
    pub fn new(extra_percent: usize) -> Self {
        Self {
            extra_percent,
            indexed: IndexedCandidates::default(),
            positions: Positions::default(),
            protected: FxHashSet::default(),
            order: Vec::new(),
            scratch: ActionScratch::new(),
            generation_batch: 256,
            outputs: Vec::new(),
        }
    }
}
impl<const C: usize> DeferredPruningProposer<C> {
    /// Bound parallel generation by `batch_sources`; one selects serial generation.
    pub fn with_generation_batch(extra_percent: usize, batch_sources: usize) -> Self {
        Self {
            generation_batch: batch_sources.max(1),
            ..Self::new(extra_percent)
        }
    }
}
impl<const C: usize> Default for DeferredPruningProposer<C> {
    fn default() -> Self {
        Self::new(25)
    }
}
impl<const C: usize> CandidateProposer<C> for DeferredPruningProposer<C> {
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error> {
        self.indexed.reset_reusing_table();
        let room = maximum.saturating_sub(buffer.basis().len());
        if room == 0 {
            return Ok(());
        }
        let watermark = room.saturating_add((room.saturating_mul(self.extra_percent) / 100).max(1));
        self.positions.reset(buffer);
        self.protected.clear();
        self.protected.extend(protected.iter().copied());
        source_order(&mut self.order, buffer.coefficients());
        if self.generation_batch == 1 {
            for &index in &self.order {
                let source = buffer.basis()[index];
                let coefficient = buffer.coefficients()[index];
                self.scratch.fill(model, &source, true);
                for &(word, value) in &self.scratch.terms {
                    if !self.protected.contains(&word)
                        && self.positions.get(buffer, &word).is_none()
                    {
                        self.indexed.add(buffer, word, coefficient * value);
                    }
                }
                if self.indexed.len() > watermark {
                    self.indexed.cap(buffer, room);
                }
            }
        } else {
            self.outputs.resize_with(self.generation_batch, Vec::new);
            for sources in self.order.chunks(self.generation_batch) {
                let positions = &self.positions;
                let protected = &self.protected;
                let coefficients = buffer.coefficients();
                let view = &*buffer;
                self.outputs
                    .par_iter_mut()
                    .zip(sources.par_iter())
                    .for_each_init(ActionScratch::new, |scratch, (output, &index)| {
                        output.clear();
                        scratch.fill(model, &view.basis()[index], true);
                        for &(word, value) in &scratch.terms {
                            if !protected.contains(&word) && positions.get(view, &word).is_none() {
                                output.push((word, coefficients[index] * value));
                            }
                        }
                    });
                for output in self.outputs.iter().take(sources.len()) {
                    for &(word, value) in output {
                        self.indexed.add(buffer, word, value);
                    }
                    if self.indexed.len() > watermark {
                        self.indexed.cap(buffer, room);
                    }
                }
            }
            for output in &mut self.outputs {
                output.clear();
            }
        }
        self.indexed.cap(buffer, room);
        Ok(())
    }
    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
        let count = self.indexed.select(buffer, threshold);
        self.indexed.reset_reusing_table();
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{JumpInput, parse_pauli_string, word_from_codes};
    use fxhash::FxHashMap;
    use num::Complex;

    fn selected<const C: usize, P: CandidateProposer<C>>(
        model: &LindbladSpec<C>,
        words: &[Word<C>],
        coefficients: &[f64],
        protected: &[Word<C>],
        room: usize,
        threshold: f64,
        proposer: &mut P,
    ) -> Vec<Word<C>> {
        let mut words = words.to_vec();
        let mut numeric = Vec::new();
        let mut view = CandidateBuffer::new(&mut words, coefficients, &mut numeric);
        proposer
            .propose(model, &mut view, protected, coefficients.len() + room)
            .unwrap();
        let count = proposer.select(&mut view, threshold);
        let result = view.candidates()[..count].to_vec();
        assert!(count <= room);
        assert_eq!(
            result.iter().copied().collect::<FxHashSet<_>>().len(),
            count
        );
        assert!(
            result
                .iter()
                .all(|word| !view.basis().contains(word) && !protected.contains(word))
        );
        result
    }

    #[test]
    fn greedy_matches_sorted_cartesian_product_without_merging() {
        let h = vec![
            ("XII".into(), 0.31),
            ("IYY".into(), -0.73),
            ("ZXI".into(), 1.17),
        ];
        let model = LindbladSpec::<2>::new(3, &h, &[]).unwrap();
        let words: Vec<_> = [1, 5, 9, 14, 31, 46]
            .iter()
            .map(|&code| {
                word_from_codes(
                    &(0..3)
                        .map(|q| ((code >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let coefficients = [0.2, -0.9, 0.13, 0.7, -0.11, 0.4];
        let mut terms: Vec<_> = model.hamiltonian_terms().collect();
        terms.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
        let mut order = Vec::new();
        source_order(&mut order, &coefficients);
        let mut pairs = Vec::new();
        for (source_rank, &source) in order.iter().enumerate() {
            for (h_rank, &(term, coefficient)) in terms.iter().enumerate() {
                let (word, sign) = comm_product(&term, &words[source]);
                if sign != 0.0 && !words.contains(&word) {
                    pairs.push((
                        word,
                        (coefficient * coefficients[source]).abs(),
                        source_rank,
                        h_rank,
                    ));
                }
            }
        }
        pairs.sort_by(|a, b| {
            b.1.total_cmp(&a.1)
                .then_with(|| a.2.cmp(&b.2))
                .then_with(|| a.3.cmp(&b.3))
        });
        let mut expected = Vec::new();
        for (word, _, _, _) in pairs {
            if !expected.contains(&word) {
                expected.push(word);
            }
        }
        let mut greedy = GreedyProposer::default();
        for room in [0, 1, 3, 100] {
            let result = selected(&model, &words, &coefficients, &[], room, 0.0, &mut greedy);
            let mut actual = result.clone();
            actual.sort_unstable();
            let mut wanted = expected[..room.min(expected.len())].to_vec();
            wanted.sort_unstable();
            assert_eq!(actual, wanted);
            assert_eq!(
                result,
                selected(
                    &model,
                    &words,
                    &coefficients,
                    &[],
                    room,
                    0.0,
                    &mut BatchedProposer::new(1)
                )
            );
        }
    }

    #[test]
    fn batch_budget_carries_unspent_slots_forward() {
        let model = LindbladSpec::<2>::new(3, &[("ZII".into(), 1.0)], &[]).unwrap();
        let words: Vec<_> = ["IZI", "IIZ", "XZI", "XIZ"]
            .iter()
            .map(|s| parse_pauli_string(s, 3).unwrap().0)
            .collect();
        let result = selected(
            &model,
            &words,
            &[8.0, 7.0, 2.0, 1.0],
            &[],
            2,
            0.0,
            &mut BatchedProposer::new(2),
        );
        assert_eq!(
            result,
            ["YZI", "YIZ"]
                .iter()
                .map(|s| parse_pauli_string(s, 3).unwrap().0)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn bounded_greedy_matches_exhaustive_max_contributions_across_replacements_and_ties() {
        let mut proposer = GreedyProposer::<2>::default();
        for trial in 0..20 {
            let h: Vec<_> = (0..12)
                .map(|i| {
                    let code = (i * 73 + trial * 11 + 1) % 256;
                    let word = (0..4)
                        .map(|q| ['I', 'X', 'Z', 'Y'][((code >> (2 * q)) & 3) as usize])
                        .collect();
                    (
                        word,
                        if trial % 2 == 0 {
                            1.0
                        } else {
                            0.2 + (i as f64 * 1.7).sin()
                        },
                    )
                })
                .collect();
            let model = LindbladSpec::new(4, &h, &[]).unwrap();
            let words: Vec<_> = (0..32)
                .map(|i| {
                    word_from_codes(
                        &(0..4)
                            .map(|q| (((i * 37 + trial * 7) >> (2 * q)) & 3) as u8)
                            .collect::<Vec<_>>(),
                    )
                    .unwrap()
                })
                .collect();
            let coefficients: Vec<_> = (0..32)
                .map(|i| {
                    if trial % 2 == 0 {
                        if i % 3 == 0 { -1.0 } else { 1.0 }
                    } else {
                        (i as f64 * 1.13).cos()
                    }
                })
                .collect();
            let mut maxima = FxHashMap::default();
            for (word, &coefficient) in words.iter().zip(&coefficients) {
                for (h_word, h_coefficient) in model.hamiltonian_terms() {
                    let (output, sign) = comm_product(&h_word, word);
                    if sign != 0.0 && !words.contains(&output) {
                        let score = (sign * h_coefficient * coefficient).abs();
                        let entry = maxima.entry(output).or_insert(0.0_f64);
                        *entry = entry.max(score);
                    }
                }
            }
            let mut scores: Vec<_> = maxima.values().copied().collect();
            scores.sort_unstable_by(|a, b| b.total_cmp(a));
            for room in [1, 3, 17, 100] {
                let result = selected(&model, &words, &coefficients, &[], room, 0.0, &mut proposer);
                let mut actual: Vec<_> = result.iter().map(|word| maxima[word]).collect();
                actual.sort_unstable_by(|a, b| b.total_cmp(a));
                assert_eq!(actual.len(), room.min(scores.len()));
                for (a, b) in actual.iter().zip(&scores) {
                    assert!((a - b).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn deferred_parallel_batches_preserve_serial_selection() {
        let model = LindbladSpec::<2>::new(
            5,
            &[
                ("XXIII".into(), 0.3),
                ("IYYII".into(), -0.7),
                ("IIXXI".into(), 0.6),
                ("IIIYY".into(), 0.5),
            ],
            &[],
        )
        .unwrap();
        let words: Vec<_> = (0..512)
            .map(|i| {
                word_from_codes(
                    &(0..5)
                        .map(|q| (((i * 37) >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let coefficients: Vec<_> = (0..512).map(|i| (i as f64 * 1.17).sin()).collect();
        for room in [8, 32, 100] {
            let serial = selected(
                &model,
                &words,
                &coefficients,
                &[],
                room,
                0.0,
                &mut DeferredPruningProposer::with_generation_batch(25, 1),
            );
            for batch in [32, 64, 256] {
                assert_eq!(
                    serial,
                    selected(
                        &model,
                        &words,
                        &coefficients,
                        &[],
                        room,
                        0.0,
                        &mut DeferredPruningProposer::with_generation_batch(25, batch)
                    )
                );
            }
        }
    }

    #[test]
    fn exact_shortlist_scoring_detects_complete_cancellation() {
        let model =
            LindbladSpec::<2>::new(2, &[("XI".into(), 1.0), ("IX".into(), 1.0)], &[]).unwrap();
        let words: Vec<_> = ["ZY", "YZ"]
            .iter()
            .map(|s| parse_pauli_string(s, 2).unwrap().0)
            .collect();
        assert_eq!(
            selected(
                &model,
                &words,
                &[1.0, -1.0],
                &[],
                2,
                0.0,
                &mut GreedyProposer::default()
            )
            .len(),
            2
        );
        assert!(
            selected(
                &model,
                &words,
                &[1.0, -1.0],
                &[],
                2,
                0.0,
                &mut RefinedGreedyProposer::default()
            )
            .is_empty()
        );
        assert!(
            selected(
                &model,
                &words,
                &[1.0, -1.0],
                &[],
                2,
                0.0,
                &mut DeferredPruningProposer::default()
            )
            .is_empty()
        );
    }

    #[test]
    fn complete_shortlist_matches_full_signed_leakage_with_duplicate_h_terms() {
        let model = LindbladSpec::<2>::new(
            4,
            &[
                ("XXII".into(), 0.3),
                ("IYYI".into(), -0.7),
                ("ZIIX".into(), 0.2),
                ("XXII".into(), -0.1),
            ],
            &[],
        )
        .unwrap();
        let words: Vec<_> = (0..50)
            .map(|i| {
                word_from_codes(
                    &(0..4)
                        .map(|q| (((i * 37) >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let coefficients: Vec<_> = (0..50).map(|i| (i as f64 * 1.71).sin()).collect();
        let mut exact = FxHashMap::default();
        for (word, &coefficient) in words.iter().zip(&coefficients) {
            for (output, value) in model.action(word) {
                if !words.contains(&output) {
                    *exact.entry(output).or_insert(0.0_f64) += coefficient * value;
                }
            }
        }
        let mut values: Vec<_> = exact.values().map(|value| value * value).collect();
        values.sort_by(|a, b| b.total_cmp(a));
        for threads in [1, 2, 6] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let result = pool.install(|| {
                selected(
                    &model,
                    &words,
                    &coefficients,
                    &[],
                    5,
                    0.0,
                    &mut RefinedGreedyProposer::new(100_000),
                )
            });
            let power: f64 = result.iter().map(|word| exact[word].powi(2)).sum();
            assert!((power - values[..5].iter().sum::<f64>()).abs() < 1e-12);
        }
    }

    #[test]
    fn hamiltonian_proposers_preserve_general_jump_leakage_by_fallback() {
        let jump = JumpInput {
            lincomb: vec![
                ("XI".into(), Complex::new(0.5, 0.0)),
                ("YI".into(), Complex::new(0.0, 0.5)),
            ],
            rate: 0.2,
        };
        let model = LindbladSpec::<2>::new(2, &[("XX".into(), 0.3)], &[jump]).unwrap();
        let words = vec![parse_pauli_string("ZI", 2).unwrap().0];
        let expected = selected(
            &model,
            &words,
            &[1.0],
            &[],
            8,
            0.0,
            &mut StreamingProposer::default(),
        );
        assert!(!expected.is_empty());
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut GreedyProposer::default()
            )
        );
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut BatchedProposer::default()
            )
        );
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut RefinedGreedyProposer::default()
            )
        );
    }

    #[test]
    fn strategies_support_wide_words_and_protected_external_strings() {
        let mut x = vec![0; 129];
        x[128] = 1;
        x[0] = 1;
        let mut z = vec![0; 129];
        z[128] = 2;
        let h = x
            .iter()
            .map(|code| ['I', 'X', 'Z', 'Y'][*code as usize])
            .collect::<String>();
        let model = LindbladSpec::<4>::new(129, &[(h, 0.3)], &[]).unwrap();
        let words = vec![word_from_codes(&z).unwrap()];
        let expected = selected(
            &model,
            &words,
            &[1.0],
            &[],
            8,
            0.0,
            &mut StreamingProposer::default(),
        );
        assert_eq!(expected.len(), 1);
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut GreedyProposer::default()
            )
        );
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut BatchedProposer::default()
            )
        );
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut RefinedGreedyProposer::default()
            )
        );
        assert_eq!(
            expected,
            selected(
                &model,
                &words,
                &[1.0],
                &[],
                8,
                0.0,
                &mut DeferredPruningProposer::default()
            )
        );
        assert!(
            selected(
                &model,
                &words,
                &[1.0],
                &expected,
                8,
                0.0,
                &mut GreedyProposer::default()
            )
            .is_empty()
        );
    }
}
