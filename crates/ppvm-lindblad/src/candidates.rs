// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! First- and second-hop candidate proposal for the real predictor-corrector
//! step, accumulated in the spare tail of the main basis array.
//!
//! [`StreamingProposer`] computes the same candidates, rates and admission
//! order as [`LindbladSpec::leakage_with_prune`] followed by admission, but
//! stores candidate words in the basis array's tail and their rates in a
//! reused numeric buffer instead of a separate word-keyed hash map.

use crate::direct_expm::ActionScratch;
use crate::indexed_candidates::IndexedCandidates;
use crate::truncate::order_by_desc_mag;
use crate::word::{Word, word_hash};
use crate::{Error, LindbladSpec};
use fxhash::FxHashMap;
use rayon::prelude::*;

/// Sources per room cap, as in [`LindbladSpec::leakage_with_prune`].
const CHUNK_SIZE: usize = 4096;
/// Sources whose actions are generated in parallel before serial accumulation.
/// Bounds the temporary output words independently of `CHUNK_SIZE`.
const GENERATION_BATCH: usize = 2048;

/// The active basis plus its spare tail, where candidates are accumulated.
///
/// `words[..prefix]` is the active basis and stays unchanged during proposal;
/// `words[prefix..]` holds candidates. `numeric` holds one rate per candidate
/// slot. Dropping the buffer without [`StreamingProposer::admit`] discards the
/// tail.
pub(crate) struct CandidateBuffer<'a, const C: usize> {
    pub(crate) words: &'a mut Vec<Word<C>>,
    pub(crate) numeric: &'a mut Vec<f64>,
    pub(crate) prefix: usize,
    coefficients: &'a [f64],
}

impl<'a, const C: usize> CandidateBuffer<'a, C> {
    /// View `words` as the active basis with coefficients `coefficients`.
    /// `numeric` is cleared and reused for candidate rates.
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

    pub(crate) fn basis(&self) -> &[Word<C>] {
        &self.words[..self.prefix]
    }

    pub(crate) fn coefficients(&self) -> &[f64] {
        self.coefficients
    }
}

impl<const C: usize> Drop for CandidateBuffer<'_, C> {
    fn drop(&mut self) {
        self.words.truncate(self.prefix);
        self.numeric.clear();
    }
}

/// Chunk-capped signed leakage accumulation.
///
/// Sources are processed in descending coefficient magnitude, and after every
/// `CHUNK_SIZE` sources only the `room` largest candidates are kept. Signed
/// sums and hash-table visitation order match
/// [`LindbladSpec::leakage_with_prune`], so admission and ties are identical.
#[derive(Default)]
pub(crate) struct StreamingProposer<const C: usize> {
    indexed: IndexedCandidates,
    outputs: Vec<Vec<(Word<C>, f64)>>,
}

impl<const C: usize> StreamingProposer<C> {
    /// Accumulate the leakage of `buffer`'s active basis into its tail.
    /// `maximum` bounds the enlarged basis, including the active prefix.
    pub(crate) fn propose(
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
                            scratch.fill(model, &basis[index], true);
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

    /// Append the candidates whose rate magnitude exceeds `threshold` to the
    /// active basis. The caller extends the coefficients with zeros.
    pub(crate) fn admit(&mut self, mut buffer: CandidateBuffer<'_, C>, threshold: f64) {
        let count = self.indexed.select(&mut buffer, threshold);
        self.indexed.reset();
        buffer.prefix += count;
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
        buffer.prefix += count;
        drop(buffer);
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
        proposer.admit(buffer, f64::NAN);
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
                    proposer.admit(buffer, 0.0);
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
