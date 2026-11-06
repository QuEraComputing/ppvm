// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Candidate proposal using a borrowed view of the adaptive state's spare storage.

use crate::indexed_candidates::IndexedCandidates;
use crate::truncate::order_by_desc_mag;
use crate::word::{Word, word_hash};
use crate::{Error, LindbladSpec};
use fxhash::FxHashMap;
use rayon::prelude::*;

const CHUNK_SIZE: usize = 4096;
const GENERATION_BATCH: usize = 2048;

/// Borrowed candidate storage with a read-only active prefix.
///
/// Candidate words occupy the main word array's tail. Numeric scratch is reused
/// by the exponential after admission. Dropping an unfinished view discards its
/// tail without changing the active words or coefficients.
pub struct CandidateBuffer<'a, const C: usize = { crate::W_CHUNKS }> {
    pub(crate) words: &'a mut Vec<Word<C>>,
    pub(crate) numeric: &'a mut Vec<f64>,
    pub(crate) prefix: usize,
    source: Option<&'a [Word<C>]>,
    coefficients: &'a [f64],
}

impl<'a, const C: usize> CandidateBuffer<'a, C> {
    /// View the active sum and use its spare word capacity for candidates.
    /// Numeric scratch is cleared while retaining its allocation.
    ///
    /// # Panics
    /// Panics when the active word and coefficient lengths differ.
    pub fn new(
        words: &'a mut Vec<Word<C>>,
        coefficients: &'a [f64],
        numeric: &'a mut Vec<f64>,
    ) -> Self {
        assert_eq!(words.len(), coefficients.len());
        let prefix = words.len();
        numeric.clear();
        Self {
            words,
            numeric,
            prefix,
            source: None,
            coefficients,
        }
    }

    pub(crate) fn separate(
        source: &'a [Word<C>],
        coefficients: &'a [f64],
        words: &'a mut Vec<Word<C>>,
        numeric: &'a mut Vec<f64>,
    ) -> Self {
        assert_eq!(source.len(), coefficients.len());
        words.clear();
        numeric.clear();
        Self {
            words,
            numeric,
            prefix: 0,
            source: Some(source),
            coefficients,
        }
    }

    /// Active strings, whose order and values remain fixed during proposal.
    pub fn basis(&self) -> &[Word<C>] {
        self.source.unwrap_or(&self.words[..self.prefix])
    }

    /// Coefficients corresponding to [`Self::basis`].
    pub fn coefficients(&self) -> &[f64] {
        self.coefficients
    }

    /// Candidate words in the writable tail, excluding the active prefix.
    pub fn candidates(&self) -> &[Word<C>] {
        &self.words[self.prefix..]
    }

    /// Mutable candidate words; the active prefix is not exposed for mutation.
    pub fn candidates_mut(&mut self) -> &mut [Word<C>] {
        &mut self.words[self.prefix..]
    }

    /// Append a candidate word, growing the backing array only when necessary.
    /// Rates and any ranking metadata are managed by the proposer.
    pub fn push(&mut self, word: Word<C>) {
        if self.words.len() == self.words.capacity() {
            self.words
                .reserve_exact((self.words.capacity() / 4).max(32));
        }
        self.words.push(word);
    }

    /// Reusable numeric storage for rates, ranking and permutation scratch.
    /// Its contents have no meaning after the proposal has been admitted.
    pub fn numeric_scratch(&mut self) -> &mut Vec<f64> {
        self.numeric
    }

    pub(crate) fn scoring_parts(&mut self) -> (&[Word<C>], &[Word<C>], &[f64], &mut Vec<f64>) {
        let (basis, candidates) = self.words.split_at(self.prefix);
        (
            self.source.unwrap_or(basis),
            candidates,
            self.coefficients,
            self.numeric,
        )
    }

    pub(crate) fn finish(mut self, count: usize) {
        assert!(count <= self.candidates().len());
        self.prefix += count;
        self.words.truncate(self.prefix);
    }
}

impl<const C: usize> Drop for CandidateBuffer<'_, C> {
    fn drop(&mut self) {
        self.words.truncate(self.prefix);
        self.numeric.clear();
    }
}

/// Strategy for proposing external strings for a real predictor-corrector step.
///
/// Selected words must be distinct and outside the active basis and protected
/// set, with at most `maximum.saturating_sub(buffer.basis().len())` survivors. Implementations
/// may retain their own metadata between calls, but cannot retain the borrowed
/// buffer. Proposal and selection are separate to preserve per-phase timings.
/// Strategies must be sendable because a step may execute in a configured thread pool.
pub trait CandidateProposer<const C: usize = { crate::W_CHUNKS }>: Send {
    /// Accumulate candidate information using the supplied tail and scratch.
    /// `maximum` bounds the total enlarged basis, including the active prefix.
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error>;

    /// Compact selected words into the beginning of the candidate tail and return
    /// their count. Only rates with magnitude strictly above `threshold` qualify.
    /// Admission initializes the corresponding coefficients to zero.
    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize;
}

/// Chunk-capped leakage proposal with position-indexed candidate accumulation.
///
/// Sources are processed in descending coefficient magnitude. Signed sums and
/// hashtable visitation preserve deterministic admission and tie ordering.
#[derive(Default)]
pub struct StreamingProposer<const C: usize = { crate::W_CHUNKS }> {
    pub(crate) indexed: IndexedCandidates,
    pub(crate) outputs: Vec<Vec<(Word<C>, f64)>>,
}

impl<const C: usize> CandidateProposer<C> for StreamingProposer<C> {
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error> {
        self.indexed.reset();
        let in_basis: FxHashMap<u64, ()> =
            buffer.basis().iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();
        let order = order_by_desc_mag(buffer.coefficients());
        let room = maximum.saturating_sub(buffer.basis().len());
        self.outputs
            .resize_with(GENERATION_BATCH.div_ceil(64), Vec::new);
        for chunk in order.chunks(CHUNK_SIZE) {
            // Batches bound temporary words without changing the source chunk's pruning cadence.
            for sources in chunk.chunks(GENERATION_BATCH) {
                let basis = buffer.basis();
                let coefficients = buffer.coefficients();
                self.outputs
                    .par_iter_mut()
                    .zip(sources.par_chunks(64))
                    .for_each_init(
                        crate::direct_expm::ActionScratch::new,
                        |scratch, (output, indices)| {
                            output.clear();
                            for &index in indices {
                                scratch.fill(model, &basis[index], true);
                                for &(hash, word, value) in &scratch.terms {
                                    if !in_basis.contains_key(&hash)
                                        && !protected_set.contains_key(&hash)
                                    {
                                        output.push((word, coefficients[index] * value));
                                    }
                                }
                            }
                        },
                    );
                // Serial source order preserves coherent sums and admission ties.
                for output in self.outputs.iter().take(sources.len().div_ceil(64)) {
                    for &(word, value) in output {
                        self.indexed.add(buffer, word, value);
                    }
                }
            }
            self.indexed.cap(buffer, room);
        }
        for output in &mut self.outputs {
            output.clear();
        }
        Ok(())
    }

    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize {
        let count = self.indexed.select(buffer, threshold);
        self.indexed.reset();
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::truncate::cap_map_to_room;
    use crate::word_from_codes;

    #[test]
    fn borrowed_candidates_preserve_reference_caps_ties_and_visitation() {
        let seed = word_from_codes::<2>(&[1; 6]).unwrap();
        let mut words = vec![seed];
        let mut numeric = Vec::new();
        let mut buffer = CandidateBuffer::new(&mut words, &[0.5], &mut numeric);
        let mut indexed = IndexedCandidates::default();
        let mut reference = FxHashMap::default();
        for chunk in 0..12 {
            for i in 0..500 {
                let code = (i * 37 + chunk * 11) % 4096;
                let word = word_from_codes::<2>(
                    &(0..6)
                        .map(|q| ((code >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap();
                let value = ((i + chunk) % 9) as f64 - 4.0;
                *reference.entry(word).or_insert(0.0) += value;
                indexed.add(&mut buffer, word, value);
            }
            cap_map_to_room(&mut reference, 300);
            indexed.cap(&mut buffer, 300);
            assert_eq!(buffer.basis(), &[seed]);
        }
        let expected: Vec<_> = reference
            .into_iter()
            .filter(|(_, value)| value.abs() > 0.0)
            .collect();
        assert_eq!(
            indexed.selected_rates(&buffer, 0.0).collect::<Vec<_>>(),
            expected
        );
        let count = indexed.select(&mut buffer, 1.5);
        buffer.finish(count);
        assert_eq!(
            words,
            std::iter::once(seed)
                .chain(
                    expected
                        .iter()
                        .filter(|(_, value)| value.abs() > 1.5)
                        .map(|(word, _)| *word)
                )
                .collect::<Vec<_>>()
        );
        assert!(numeric.is_empty());
    }

    #[test]
    fn nonfinite_threshold_rejects_tail_and_preserves_prefix() {
        let seed = word_from_codes::<2>(&[1, 0]).unwrap();
        let candidate = word_from_codes::<2>(&[0, 2]).unwrap();
        let mut words = vec![seed];
        let mut numeric = Vec::new();
        let mut buffer = CandidateBuffer::new(&mut words, &[1.0], &mut numeric);
        let mut indexed = IndexedCandidates::default();
        indexed.add(&mut buffer, candidate, 0.5);
        let count = indexed.select(&mut buffer, f64::NAN);
        buffer.finish(count);
        assert_eq!(words, vec![seed]);
    }
}
