// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::SmallRng;
#[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
use rayon::prelude::*;

use super::GeneralizedTableauMixture;
use crate::{Bitstring, GeneralizedTableau, MeasureScratch, RowStorage};

/// Immutable mixture snapshot plus an independent seeded shot stream.
#[derive(Clone)]
pub struct MixtureSampler<A: RowStorage, I: Bitstring, H> {
    pub entries: Vec<(GeneralizedTableau<A, I, H>, f64)>,
    cumulative: Vec<f64>,
    rng: SmallRng,
    scratch: MeasureScratch<I>,
}

impl<A: RowStorage, I: Bitstring, H> MixtureSampler<A, I, H> {
    pub(crate) fn new(
        entries: Vec<(GeneralizedTableau<A, I, H>, f64)>,
        cumulative: Vec<f64>,
        rng: SmallRng,
    ) -> Self {
        Self {
            entries,
            cumulative,
            rng,
            scratch: MeasureScratch::new(),
        }
    }

    fn choice(&mut self) -> (usize, u64) {
        let probability = self.rng.random::<f64>();
        let index = self
            .cumulative
            .partition_point(|&bound| bound <= probability)
            .min(self.entries.len().saturating_sub(1));
        (index, self.rng.random())
    }

    pub fn sample(&mut self) -> Vec<Option<bool>> {
        let (index, seed) = self.choice();
        let mut tab = self.entries[index].0.fork(Some(seed));
        tab.measure_all_with_scratch(&mut self.scratch)
    }

    pub fn sample_shots_serial(&mut self, shots: usize) -> Vec<Vec<Option<bool>>> {
        (0..shots).map(|_| self.sample()).collect()
    }

    #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
    pub fn sample_shots_parallel(&mut self, shots: usize) -> Vec<Vec<Option<bool>>>
    where
        A: Send + Sync,
        I: Send + Sync,
        H: Sync,
    {
        let choices: Vec<_> = (0..shots).map(|_| self.choice()).collect();
        choices
            .par_iter()
            .map_init(MeasureScratch::new, |scratch, &(index, seed)| {
                let mut tab = self.entries[index].0.fork(Some(seed));
                tab.measure_all_with_scratch(scratch)
            })
            .collect()
    }

    #[cfg(any(not(feature = "rayon"), target_arch = "wasm32"))]
    pub fn sample_shots(&mut self, shots: usize) -> Vec<Vec<Option<bool>>> {
        self.sample_shots_serial(shots)
    }

    #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
    pub fn sample_shots(&mut self, shots: usize) -> Vec<Vec<Option<bool>>>
    where
        A: Send + Sync,
        I: Send + Sync,
        H: Sync,
    {
        let threads = rayon::current_num_threads();
        if threads <= 1 || shots < 4 * threads {
            self.sample_shots_serial(shots)
        } else {
            self.sample_shots_parallel(shots)
        }
    }
}

impl<A, I, H> GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    pub fn sampler(&mut self) -> MixtureSampler<A, I, H> {
        let mut entries = self.entries.clone();
        entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut sum = 0.0;
        let cumulative = entries
            .iter()
            .map(|entry| {
                sum += entry.1;
                sum
            })
            .collect();
        let rng = SmallRng::seed_from_u64(self.rng.random());
        MixtureSampler::new(entries, cumulative, rng)
    }
}
