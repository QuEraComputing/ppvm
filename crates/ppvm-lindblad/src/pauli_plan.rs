// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared parity tables and local phase rules for Pauli Hamiltonian actions.

use crate::algebra::comm_product;
use crate::word::{CHUNK_BITS, Word, word_from_codes, word_support};
use fxhash::FxHashSet;
use ppvm_traits::PauliWordTrait;

struct Term<const C: usize> {
    word: Word<C>,
    coefficient: f64,
    sites: Option<[usize; 2]>,
    signs: [i8; 16],
}

pub(crate) struct PauliActionPlan<const C: usize> {
    terms: Vec<Term<C>>,
    jumps: Vec<f64>,
    words: usize,
    tables: Vec<u64>,
}

fn label<const C: usize>(word: &Word<C>, site: usize) -> usize {
    let chunk = site / CHUNK_BITS;
    let shift = site % CHUNK_BITS;
    ((word.xbits.data[chunk] >> shift) & 1) as usize
        | ((((word.zbits.data[chunk] >> shift) & 1) as usize) << 1)
}

impl<const C: usize> PauliActionPlan<C> {
    pub(crate) fn has_diagonal_noise(&self) -> bool {
        self.jumps.iter().any(|rate| *rate != 0.0)
    }

    pub(crate) fn uniform_scale(&self) -> Option<f64> {
        let magnitude = self
            .terms
            .iter()
            .find(|term| term.coefficient != 0.0)
            .map_or(0.5, |term| term.coefficient.abs());
        if !magnitude.is_finite()
            || !self
                .terms
                .iter()
                .all(|term| term.coefficient == 0.0 || term.coefficient.abs() == magnitude)
        {
            return None;
        }
        Some(2.0 * magnitude)
    }

    pub(crate) fn new(n: usize, h: &[(Word<C>, f64)], jumps: &[(Word<C>, f64)]) -> Option<Self> {
        let mut seen = FxHashSet::default();
        if h.iter().any(|(word, _)| !seen.insert(*word)) {
            return None;
        }
        let words = (h.len() + jumps.len()).div_ceil(64);
        let mut tables = vec![0; n * 4 * words];
        for (index, (word, _)) in h.iter().chain(jumps).enumerate() {
            for site in 0..n {
                let generator = label(word, site);
                for pauli in 1..4 {
                    if generator != 0 && generator != pauli {
                        tables[(site * 4 + pauli) * words + index / 64] |= 1u64 << (index % 64);
                    }
                }
            }
        }
        let terms = h
            .iter()
            .map(|&(word, coefficient)| {
                let mut support = Vec::new();
                word_support(&word, &mut support);
                let sites = if (1..=2).contains(&support.len()) {
                    Some([support[0] as usize, *support.last().unwrap() as usize])
                } else {
                    None
                };
                let mut signs = [0; 16];
                if let Some([a, b]) = sites {
                    for (pattern, sign) in signs.iter_mut().enumerate() {
                        let mut codes = vec![0; n];
                        codes[a] = (pattern & 3) as u8;
                        if b != a {
                            codes[b] = (pattern >> 2) as u8;
                        }
                        *sign = comm_product(&word, &word_from_codes::<C>(&codes).unwrap()).1 as i8;
                    }
                }
                Term {
                    word,
                    coefficient,
                    sites,
                    signs,
                }
            })
            .collect();
        Some(Self {
            terms,
            jumps: jumps.iter().map(|(_, rate)| *rate).collect(),
            words,
            tables,
        })
    }

    pub(crate) fn fill(
        &self,
        word: &Word<C>,
        active: &mut Vec<u64>,
        output: &mut Vec<(Word<C>, f64)>,
    ) {
        active.clear();
        active.resize(self.words, 0);
        for chunk in 0..C {
            let mut support = word.xbits.data[chunk] | word.zbits.data[chunk];
            while support != 0 {
                let site = chunk * CHUNK_BITS + support.trailing_zeros() as usize;
                support &= support - 1;
                let start = (site * 4 + label(word, site)) * self.words;
                for (bits, table) in active
                    .iter_mut()
                    .zip(&self.tables[start..start + self.words])
                {
                    *bits ^= table;
                }
            }
        }
        output.clear();
        let mut diagonal = 0.0;
        for (block, &initial) in active.iter().enumerate() {
            let mut bits = initial;
            while bits != 0 {
                let index = block * 64 + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                if index >= self.terms.len() {
                    diagonal += -2.0 * self.jumps[index - self.terms.len()];
                    continue;
                }
                let term = &self.terms[index];
                let (result, sign) = if let Some([a, b]) = term.sites {
                    let pattern = label(word, a) | if a == b { 0 } else { label(word, b) << 2 };
                    let mut result = *word;
                    for chunk in 0..C {
                        result.xbits.data[chunk] ^= term.word.xbits.data[chunk];
                        result.zbits.data[chunk] ^= term.word.zbits.data[chunk];
                    }
                    result.rehash();
                    (result, term.signs[pattern] as f64)
                } else {
                    comm_product(&term.word, word)
                };
                let value = term.coefficient * sign;
                if value != 0.0 {
                    output.push((result, value));
                }
            }
        }
        if diagonal != 0.0 {
            output.push((*word, diagonal));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{JumpInput, LindbladSpec, parse_pauli_string};
    use num::Complex;

    #[test]
    fn compiled_tables_match_every_local_pauli_and_diagonal_jump() {
        let h = vec![
            ("XXII".into(), 0.3),
            ("IYYI".into(), -0.7),
            ("ZIIX".into(), 0.2),
            ("IXYZ".into(), 0.4),
        ];
        let jump = JumpInput {
            lincomb: vec![("IZII".into(), Complex::new(1.0, 0.0))],
            rate: 0.4,
        };
        let spec = LindbladSpec::<2>::new(4, &h, &[jump]).unwrap();
        let parsed: Vec<_> = h
            .iter()
            .map(|(s, c)| (parse_pauli_string::<2>(s, 4).unwrap().0, *c))
            .collect();
        let jumps = [(parse_pauli_string::<2>("IZII", 4).unwrap().0, 0.4)];
        let plan = PauliActionPlan::new(4, &parsed, &jumps).unwrap();
        let mut active = Vec::new();
        let mut output = Vec::new();
        for code in 0..256 {
            let word = word_from_codes::<2>(
                &(0..4)
                    .map(|q| ((code >> (2 * q)) & 3) as u8)
                    .collect::<Vec<_>>(),
            )
            .unwrap();
            plan.fill(&word, &mut active, &mut output);
            let mut expected = spec.action(&word);
            output.sort_unstable_by_key(|(w, _)| *w);
            expected.sort_unstable_by_key(|(w, _)| *w);
            assert_eq!(output, expected);
        }
    }
}
