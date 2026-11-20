// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Candidate proposal for the real predictor-corrector step.
//!
//! A [`CandidateProposer`] chooses which strings outside the basis a step
//! admits. Candidates are written to the spare tail of the main basis array,
//! viewed through a [`CandidateBuffer`]. [`StreamingProposer`], the default,
//! computes the same candidates, rates and admission order as
//! [`LindbladSpec::leakage_with_prune`] followed by admission.

use crate::indexed_candidates::IndexedCandidates;
use crate::truncate::order_by_desc_mag;
use crate::word::{Word, word_hash};
use crate::{Error, LindbladSpec};
use fxhash::FxHashMap;
use num::Complex;
use rayon::prelude::*;

/// Sources per room cap, as in [`LindbladSpec::leakage_with_prune`].
const CHUNK_SIZE: usize = 4096;
/// Sources whose actions are generated in parallel before serial accumulation.
/// Bounds the temporary output words independently of `CHUNK_SIZE`.
const GENERATION_BATCH: usize = 2048;

/// Strategy that chooses which strings outside the basis a step admits.
///
/// [`LindbladSpec::pc_step`] uses [`StreamingProposer`];
/// [`LindbladSpec::pc_step_with_proposer`] takes any implementation. A step
/// calls `propose` and then `select` once per enlargement (at most twice).
pub trait CandidateProposer<const C: usize = { crate::W_CHUNKS }>: Send {
    /// Write candidate words to `buffer`'s tail with [`CandidateBuffer::push`].
    /// Candidates must be distinct and outside the basis and `protected`.
    /// `maximum` bounds the enlarged basis, including the current basis.
    fn propose(
        &mut self,
        model: &LindbladSpec<C>,
        buffer: &mut CandidateBuffer<'_, C>,
        protected: &[Word<C>],
        maximum: usize,
    ) -> Result<(), Error>;

    /// Move the admitted candidates to the start of the tail and return how
    /// many there are: at most `maximum - basis().len()`, each with a rate
    /// magnitude above `threshold`. They enter the basis with coefficient zero.
    fn select(&mut self, buffer: &mut CandidateBuffer<'_, C>, threshold: f64) -> usize;
}

/// The current basis plus its spare tail, where candidates are collected.
///
/// The basis and its coefficients stay unchanged during proposal. Candidates
/// the step does not admit are discarded when the buffer is dropped.
pub struct CandidateBuffer<'a, const C: usize = { crate::W_CHUNKS }> {
    pub(crate) words: &'a mut Vec<Word<C>>,
    pub(crate) numeric: &'a mut Vec<f64>,
    pub(crate) prefix: usize,
    coefficients: &'a [f64],
}

impl<'a, const C: usize> CandidateBuffer<'a, C> {
    /// View `words` as the basis with coefficients `coefficients`. `numeric`
    /// is cleared and lent to the proposer as scratch.
    pub(crate) fn new(
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
            coefficients,
        }
    }

    /// The current basis.
    pub fn basis(&self) -> &[Word<C>] {
        &self.words[..self.prefix]
    }

    /// Coefficients of [`Self::basis`].
    pub fn coefficients(&self) -> &[f64] {
        self.coefficients
    }

    /// Candidate words written so far.
    pub fn candidates(&self) -> &[Word<C>] {
        &self.words[self.prefix..]
    }

    /// Candidate words, for reordering in [`CandidateProposer::select`].
    pub fn candidates_mut(&mut self) -> &mut [Word<C>] {
        &mut self.words[self.prefix..]
    }

    /// Append a candidate word. The tail shares the main basis array, so it
    /// grows by a quarter of its capacity instead of doubling it.
    pub fn push(&mut self, word: Word<C>) {
        if self.words.len() == self.words.capacity() {
            self.words
                .reserve_exact((self.words.capacity() / 4).max(32));
        }
        self.words.push(word);
    }

    /// Scratch the proposer may use for rates and ranking. The step reuses
    /// it afterwards, so its contents do not survive `select`.
    pub fn numeric_scratch(&mut self) -> &mut Vec<f64> {
        self.numeric
    }

    /// Admit the first `count` candidates into the basis.
    pub(crate) fn finish(mut self, count: usize) {
        assert!(count <= self.candidates().len());
        self.prefix += count;
    }
}

impl<const C: usize> Drop for CandidateBuffer<'_, C> {
    fn drop(&mut self) {
        self.words.truncate(self.prefix);
        self.numeric.clear();
    }
}

/// Per-thread scratch for one string's generator action.
pub(crate) struct ActionScratch<const C: usize> {
    support: Vec<u32>,
    candidates: Vec<u32>,
    local: FxHashMap<Word<C>, Complex<f64>>,
    active: Vec<u64>,
    /// Output words with their [`word_hash`] and generator values, in the
    /// order of [`LindbladSpec::compute_action_terms`].
    pub(crate) terms: Vec<(u64, Word<C>, f64)>,
}

impl<const C: usize> ActionScratch<C> {
    pub(crate) fn new() -> Self {
        Self {
            support: Vec::new(),
            candidates: Vec::new(),
            local: FxHashMap::default(),
            active: Vec::new(),
            terms: Vec::new(),
        }
    }

    /// Fill [`Self::terms`] with `L*(word)`. Runs once per source string in
    /// the generation loop, whose speed depends on this being inlined.
    #[inline]
    pub(crate) fn fill(&mut self, spec: &LindbladSpec<C>, word: &Word<C>) {
        if let Some(plan) = &spec.pauli_plan {
            plan.fill_keyed(word, &mut self.active, &mut self.terms);
            self.terms
                .sort_unstable_by(|(ha, a, _), (hb, b, _)| ha.cmp(hb).then_with(|| a.cmp(b)));
            return;
        }
        let terms = spec.compute_action_terms(
            word,
            &mut self.support,
            &mut self.candidates,
            &mut self.local,
        );
        self.terms.clear();
        self.terms.extend(
            terms
                .into_iter()
                .map(|(word, value)| (word_hash(&word), word, value)),
        );
    }
}

/// The default proposer: chunk-capped signed leakage accumulation.
///
/// Sources are processed in descending coefficient magnitude, and after every
/// `CHUNK_SIZE` sources only the `room` largest candidates are kept. Signed
/// sums and hash-table visitation order match
/// [`LindbladSpec::leakage_with_prune`], so admission and ties are identical.
#[derive(Default)]
pub struct StreamingProposer<const C: usize = { crate::W_CHUNKS }> {
    indexed: IndexedCandidates,
    outputs: Vec<Vec<(Word<C>, f64)>>,
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
            for sources in chunk.chunks(GENERATION_BATCH) {
                let basis = buffer.basis();
                let coefficients = buffer.coefficients();
                self.outputs
                    .par_iter_mut()
                    .zip(sources.par_chunks(64))
                    .for_each_init(ActionScratch::new, |scratch, (output, indices)| {
                        output.clear();
                        for &index in indices {
                            scratch.fill(model, &basis[index]);
                            for &(hash, word, value) in &scratch.terms {
                                if !in_basis.contains_key(&hash)
                                    && !protected_set.contains_key(&hash)
                                {
                                    output.push((word, coefficients[index] * value));
                                }
                            }
                        }
                    });
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
    fn indexed_candidates_match_map_caps_ties_and_visitation() {
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
            .filter(|(_, value)| value.abs() > 1.5)
            .map(|(word, _)| word)
            .collect();
        let count = indexed.select(&mut buffer, 1.5);
        buffer.finish(count);
        assert_eq!(
            words,
            std::iter::once(seed).chain(expected).collect::<Vec<_>>()
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
        let mut proposer = StreamingProposer::default();
        proposer.indexed.add(&mut buffer, candidate, 0.5);
        let count = proposer.select(&mut buffer, f64::NAN);
        buffer.finish(count);
        assert_eq!(words, vec![seed]);
    }

    /// Several room caps, more than one `CHUNK_SIZE` of sources, and thread
    /// counts that split generation batches differently.
    #[test]
    fn streaming_admission_matches_leakage_with_prune() {
        let n = 7;
        let mut h = Vec::new();
        for site in 0..n - 1 {
            for symbol in ['X', 'Y'] {
                let mut term = vec!['I'; n];
                term[site] = symbol;
                term[site + 1] = symbol;
                h.push((term.into_iter().collect(), 0.3));
            }
        }
        let spec = LindbladSpec::<2>::new(n, &h, &[]).unwrap();
        let basis: Vec<_> = (0..14000)
            .filter(|i| i % 3 != 0)
            .map(|i| {
                word_from_codes::<2>(
                    &(0..n)
                        .map(|q| ((i >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let coefficients: Vec<_> = (0..basis.len()).map(|i| (i as f64).sin()).collect();
        for room in [0, 100, 2000] {
            let maximum = basis.len() + room;
            let expected: Vec<_> = spec
                .leakage_with_prune(&basis, &coefficients, &[], maximum, 0.0)
                .unwrap()
                .into_iter()
                .map(|(word, _)| word)
                .collect();
            for threads in [1, 2, 6] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let mut words = basis.clone();
                let mut numeric = Vec::new();
                pool.install(|| {
                    let mut buffer = CandidateBuffer::new(&mut words, &coefficients, &mut numeric);
                    let mut proposer = StreamingProposer::default();
                    proposer.propose(&spec, &mut buffer, &[], maximum).unwrap();
                    let count = proposer.select(&mut buffer, 0.0);
                    buffer.finish(count);
                });
                assert_eq!(
                    words[basis.len()..],
                    expected[..],
                    "room={room}, threads={threads}"
                );
            }
        }
    }
}
