// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::marker::PhantomData;

use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};
use ppvm_tableau::{
    data::GeneralizedTableau, prelude::Config, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use crate::{
    sampler::Sampler,
    storage::{EntryStore, phase_loss_hash, vec::VecStorage, word_fingerprint},
};

#[derive(Clone)]
pub struct GeneralizedTableauSum<
    T: Config,
    I,
    C: SparseVector<Complex<T::Coeff>, I> = Vec<(Complex64, I)>,
    S: EntryStore<T, I, C> = VecStorage<T, I, C>,
> {
    pub n_qubits: usize,
    pub entries: S,
    pub(crate) rng: SmallRng,
    pub(crate) sum_cutoff: T::Coeff,
    _phantom: PhantomData<(I, C)>,
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>>
    GeneralizedTableauSum<T, I, C, S>
where
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
{
    pub fn new(n_qubits: usize, coefficient_threshold: T::Coeff, sum_cutoff: T::Coeff) -> Self {
        let rng = rand::make_rng();
        let g_tab: GeneralizedTableau<T, I, C> =
            GeneralizedTableau::new(n_qubits, coefficient_threshold);
        let wfp = word_fingerprint(&g_tab);
        let plh = phase_loss_hash(&g_tab);
        let mut storage = S::with_capacity(1);
        storage.insert_or_merge_batch(vec![(g_tab, T::Coeff::one(), wfp, plh)], &sum_cutoff);
        Self {
            n_qubits,
            entries: storage,
            rng,
            sum_cutoff,
            _phantom: PhantomData,
        }
    }

    pub fn new_with_seed(
        n_qubits: usize,
        coefficient_threshold: T::Coeff,
        sum_cutoff: T::Coeff,
        seed: u64,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let tab_seed = rng.random::<u64>();
        let g_tab: GeneralizedTableau<T, I, C> =
            GeneralizedTableau::new_with_seed(n_qubits, coefficient_threshold, tab_seed);
        let wfp = word_fingerprint(&g_tab);
        let plh = phase_loss_hash(&g_tab);
        let mut storage = S::with_capacity(1);
        storage.insert_or_merge_batch(vec![(g_tab, T::Coeff::one(), wfp, plh)], &sum_cutoff);
        Self {
            n_qubits,
            entries: storage,
            rng,
            sum_cutoff,
            _phantom: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn truncate(&mut self) {
        let length_before_truncation = self.entries.len();
        let cutoff = self.sum_cutoff.clone();
        // Filter both vectors in lockstep so `entry_fingerprints` stays
        // aligned with `entries` after retain shifts indices.
        self.entries.retain(|_, p| *p > cutoff);
        if self.entries.len() < length_before_truncation {
            self.normalize_probabilities();
        }
    }

    pub fn normalize_probabilities(&mut self) {
        let norm = self
            .entries
            .iter()
            .fold(T::Coeff::zero(), |acc, entry| acc + entry.1.clone());
        self.entries.for_each_mut(|_, p| {
            *p = p.clone() / norm.clone();
        });
    }

    pub fn sampler(&mut self) -> Sampler<T, I, C> {
        // Sorting reorders entries; drop the parallel fingerprint cache
        // rather than permuting it (cheaper, and the sampler doesn't read
        // it). Subsequent noise ops on this sum would just recompute.
        let mut entries: Vec<(GeneralizedTableau<T, I, C>, T::Coeff)> = self
            .entries
            .iter()
            .map(|(t, c)| (t.clone(), c.clone()))
            .collect();
        entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut p_acc = T::Coeff::zero();
        let mut p_cum = Vec::<T::Coeff>::with_capacity(entries.len());
        for entry in entries.iter() {
            p_acc += entry.1.clone();
            p_cum.push(p_acc.clone())
        }

        let seed = self.rng.random::<u64>();
        let rng = SmallRng::seed_from_u64(seed);

        debug_assert!(
            *p_cum.last().unwrap_or(&T::Coeff::zero()) >= T::Coeff::one() - self.sum_cutoff.clone(),
            "Normalization error in sum"
        );

        Sampler {
            p_cumulative: p_cum,
            entries,
            rng,
            scratch: ppvm_tableau::measure::MeasureScratch::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_pauli_sum::config::fxhash::ByteF64;
    use ppvm_traits::traits::{
        Clifford, Depolarizing, LossChannel, LossyMeasure, PauliError, Reset, ResetLossChannel,
        TGate,
    };

    use crate::storage::map::MapStorage;

    type TestConfig = ByteF64<1>;
    type TestCoeffVec = Vec<(Complex64, u128)>;
    type TestSum = GeneralizedTableauSum<TestConfig, u128>;
    type TestMapSum = GeneralizedTableauSum<
        TestConfig,
        u128,
        TestCoeffVec,
        MapStorage<TestConfig, u128, TestCoeffVec>,
    >;

    fn make(n: usize) -> TestSum {
        GeneralizedTableauSum::new_with_seed(n, 1e-12, 1e-10, 42)
    }

    fn sum_of_probabilities(tab: &TestSum) -> f64 {
        tab.entries.iter().map(|e| *e.1).sum()
    }

    fn outcome_probability(tab: &TestSum, outcome: bool) -> f64 {
        tab.entries
            .iter()
            .filter_map(|(entry, p)| {
                let mut entry = entry.clone();
                (entry.measure(0) == Some(outcome)).then_some(*p)
            })
            .sum()
    }

    #[test]
    fn test_new_initial_state() {
        let tab = make(3);
        assert_eq!(tab.len(), 1);
        assert_eq!(tab.n_qubits, 3);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries.entries[0].0.is_lost.iter().all(|x| !x));
    }

    #[test]
    fn test_clifford_gates_dont_add_branches() {
        let mut tab = make(3);
        tab.h(0);
        tab.cnot(0, 1);
        tab.cz(0, 2);
        tab.s(1);
        tab.x(2);
        tab.y(0);
        tab.z(1);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_t_gate_doesnt_add_branches() {
        let mut tab = make(2);
        tab.h(0);
        tab.t(0);
        tab.t_adj(0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_probabilities() {
        let mut tab = make(2);
        let cloned = tab.entries.entries[0].0.clone();
        tab.entries.entries.push((cloned, 3.0));
        tab.entries.fingerprints.push(0);
        tab.entries.word_fingerprints.push(0);
        tab.entries.phase_loss_hashes.push(0);
        tab.entries.mark_keys_dirty();
        tab.normalize_probabilities();
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        assert!((tab.entries.entries[0].1 - 0.25).abs() < 1e-12);
        assert!((tab.entries.entries[1].1 - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_truncate_removes_low_probability_entries() {
        // sum_cutoff = 0.1
        let mut tab: TestSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.1, 42);
        let cloned = tab.entries.entries[0].0.clone();
        tab.entries.entries.push((cloned, 0.05));
        tab.entries.fingerprints.push(0);
        tab.entries.word_fingerprints.push(0);
        tab.entries.phase_loss_hashes.push(0);
        tab.entries.mark_keys_dirty();
        tab.normalize_probabilities();
        // entries are now ~ (0.952, 0.048); 0.048 is below the 0.1 cutoff
        tab.truncate();
        assert_eq!(tab.len(), 1);
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_loss_channel_zero_probability_doesnt_branch() {
        let mut tab = make(2);
        tab.loss_channel(0, 0.0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(!tab.entries.entries[0].0.is_lost[0]);
    }

    #[test]
    fn test_loss_channel_creates_balanced_branch() {
        let mut tab = make(2);
        tab.loss_channel(0, 0.5);
        assert_eq!(tab.len(), 2);
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        let lost_count = tab.entries.iter().filter(|e| e.0.is_lost[0]).count();
        assert_eq!(lost_count, 1);
        for entry in tab.entries.iter() {
            assert!((*entry.1 - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn test_reset_case_a_branch_fingerprint_lets_branches_merge() {
        // h(0) prepares |+⟩, which makes reset(0) take the case-a path
        // (Z anticommutes with the X stabilizer). prob_0 = prob_1 = 0.5,
        // and BOTH outcome branches end up at |0⟩:
        //   - existing entry: project_case_a(outcome=false)  → |0⟩, p=0.5
        //   - new branch:     project_case_a(outcome=true), x(0) → |0⟩, p=0.5
        //
        // The two tableaux are structurally identical, so
        // `insert_or_merge_batch` MUST coalesce them into a single entry
        // with p=1.0. That coalescence relies on the pushed branch carrying
        // a phase/loss fingerprint that matches its actual post-x tableau:
        // if the fingerprint is captured BEFORE x(addr0) it hashes to the
        // wrong bucket / fp_index slot, `structurally_equal` is never
        // consulted, and the branch silently sticks around as a duplicate
        // entry (the sampling distribution stays correct, but the merge
        // optimisation is lost). This test fails (len == 2) if anyone
        // re-introduces that ordering.
        let mut tab = make(1);
        tab.h(0);
        tab.reset(0);
        assert_eq!(
            tab.len(),
            1,
            "case-a reset of |+⟩ must merge its two |0⟩ branches into one entry"
        );
        assert!(
            (tab.entries.entries[0].1 - 1.0).abs() < 1e-12,
            "merged entry must carry the full probability mass"
        );
    }

    #[test]
    fn test_reset_case_a_branch_fingerprint_lets_branches_merge_map() {
        // Same invariant as the Vec test above, but on `MapStorage`. The
        // bucket key is `word_fp ^ phase_loss`, so a stale phase/loss
        // fingerprint puts the branch in a different bucket from the
        // existing structurally-equal entry — `buckets.len()` would jump
        // to 2.
        let mut tab: TestMapSum = GeneralizedTableauSum::new_with_seed(1, 1e-12, 1e-10, 42);
        tab.h(0);
        tab.reset(0);
        assert_eq!(
            tab.len(),
            1,
            "case-a reset of |+⟩ must merge its two |0⟩ branches into one entry"
        );
        assert_eq!(
            tab.entries.buckets.len(),
            1,
            "merged entry must occupy a single fingerprint bucket"
        );
        let (_, p) = tab.entries.iter().next().unwrap();
        assert!(
            (*p - 1.0).abs() < 1e-12,
            "merged entry must carry the full probability mass"
        );
    }

    #[test]
    fn test_reset_loss_channel_merges_equal_vec_entries() {
        let mut tab = make(1);
        tab.loss_channel(0, 0.5);
        assert_eq!(tab.len(), 2);

        tab.reset_loss_channel(0);

        assert_eq!(tab.len(), 1);
        let (entry, p) = tab.entries.iter().next().unwrap();
        assert!(!entry.is_lost[0]);
        assert!((*p - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_reset_loss_channel_merges_equal_map_entries() {
        let mut tab: TestMapSum = GeneralizedTableauSum::new_with_seed(1, 1e-12, 1e-10, 42);
        tab.loss_channel(0, 0.5);
        assert_eq!(tab.len(), 2);

        tab.reset_loss_channel(0);

        assert_eq!(tab.len(), 1);
        assert_eq!(tab.entries.buckets.len(), 1);
        let (entry, p) = tab.entries.iter().next().unwrap();
        assert!(!entry.is_lost[0]);
        assert!((*p - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_depolarize_zero_probability_doesnt_branch() {
        let mut tab = make(2);
        tab.depolarize(0, 0.0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pauli_error_zero_probability_doesnt_branch() {
        let mut tab = make(1);
        tab.pauli_error(0, [0.0, 0.0, 0.0]);
        assert_eq!(tab.len(), 1);
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        assert!((outcome_probability(&tab, false) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_pauli_error_nonuniform_probabilities_on_zero_state() {
        let mut tab = make(1);
        tab.pauli_error(0, [0.2, 0.3, 0.1]);

        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        assert!((outcome_probability(&tab, true) - 0.5).abs() < 1e-12);
        assert!((outcome_probability(&tab, false) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_pauli_error_skips_already_lost_qubit() {
        let mut tab = make(1);
        tab.loss_channel(0, 1.0);
        tab.pauli_error(0, [0.2, 0.3, 0.1]);

        assert_eq!(tab.len(), 1);
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        assert!(tab.entries.entries[0].0.is_lost[0]);
    }

    #[test]
    fn test_probability_sum_preserved_under_mixed_cutoff_noise() {
        // 2 entries with imbalanced probs; depolarize creates per-entry branches.
        // The smaller entry's branches fall below sum_cutoff while the larger entry's
        // branches stay above it. The total probability mass must still sum to 1.
        let mut tab: TestSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.05, 42);
        tab.loss_channel(0, 0.9); // entries: [(orig, 0.1), (lost_q0, 0.9)]
        tab.depolarize(1, 0.3);
        let sum = sum_of_probabilities(&tab);
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "Probabilities should sum to 1, got {}",
            sum
        );
    }

    #[test]
    fn test_sampler_zero_state_always_zero() {
        let mut tab = make(3);
        let mut sampler = tab.sampler();
        for _ in 0..10 {
            let m = sampler.sample();
            assert_eq!(m.len(), 3);
            for outcome in m {
                assert_eq!(outcome, Some(false));
            }
        }
    }

    #[test]
    fn test_sampler_bell_pair_correlated() {
        let mut tab = make(2);
        tab.h(0);
        tab.cnot(0, 1);
        let mut sampler = tab.sampler();
        for _ in 0..50 {
            let m = sampler.sample();
            assert_eq!(m.len(), 2);
            assert_eq!(
                m[0], m[1],
                "Bell-pair measurements should be correlated, got {:?}",
                m
            );
        }
    }

    #[test]
    fn test_sampler_seeded_reproducible() {
        let build = || {
            let mut tab = make(2);
            tab.h(0);
            tab.cnot(0, 1);
            tab.sampler().sample_shots(20)
        };
        assert_eq!(build(), build());
    }

    #[test]
    fn test_depolarize_creates_branches_summing_to_one() {
        // depolarize on |0...0⟩ should branch into the identity + X + Y + Z branches.
        // X|0⟩ and Z|0⟩ have distinct tableau data from the original and from each other,
        // so we expect >= 2 entries and total probability == 1.
        let mut tab = make(2);
        tab.depolarize(0, 0.6);
        assert!(
            tab.len() > 1,
            "depolarize should create branches, got len={}",
            tab.len()
        );
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_loss_channel_full_loss_marks_qubit_lost() {
        // p=1.0 means certain loss: original branch's probability collapses to 0 and
        // is dropped by truncate; the lost branch is renormalized to 1.0.
        let mut tab = make(2);
        tab.loss_channel(0, 1.0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries.entries[0].0.is_lost[0]);
        assert!(!tab.entries.entries[0].0.is_lost[1]);
    }

    #[test]
    fn test_loss_channel_skips_already_lost_qubit() {
        // After a full loss on q0, applying loss to q0 again must not branch further:
        // the lone entry already has is_lost[0]=true and is skipped by the iterator.
        let mut tab = make(2);
        tab.loss_channel(0, 1.0);
        assert_eq!(tab.len(), 1);
        tab.loss_channel(0, 0.5);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries.entries[0].0.is_lost[0]);
    }

    #[test]
    fn test_depolarize_skips_already_lost_qubit() {
        // Depolarizing a lost qubit must not create branches.
        let mut tab = make(2);
        tab.loss_channel(0, 1.0);
        tab.depolarize(0, 0.3);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_sampler_returns_none_for_lost_qubit() {
        let mut tab = make(2);
        tab.loss_channel(1, 1.0);
        let mut sampler = tab.sampler();
        for _ in 0..10 {
            let m = sampler.sample();
            assert_eq!(m.len(), 2);
            assert_eq!(m[0], Some(false));
            assert_eq!(m[1], None);
        }
    }

    #[test]
    fn test_sampler_respects_sorted_branch_probabilities() {
        // loss_channel(0, 0.9) leaves the storage in insertion order
        // [(orig, 0.1), (lost, 0.9)]. The sampler sorts entries by descending
        // probability before building p_cumulative, so the cumulative array
        // must be derived from the sorted view — otherwise sampling picks
        // the wrong tableau and the lost/unlost frequencies invert.
        let mut tab = make(1);
        tab.loss_channel(0, 0.9);
        let mut sampler = tab.sampler();
        let n_shots = 4000;
        let lost_count = (0..n_shots)
            .filter(|_| sampler.sample()[0].is_none())
            .count();
        let lost_frac = lost_count as f64 / n_shots as f64;
        // Expect ~0.9; the bug produces ~0.1.
        assert!(
            (lost_frac - 0.9).abs() < 0.03,
            "expected ~90% lost shots, got {:.3}",
            lost_frac
        );
    }

    #[test]
    fn test_word_fingerprint_cache_stays_consistent() {
        // After a real gate+noise sequence (ending on a merge), the cached
        // per-entry fingerprints must stay aligned with `entries` and equal a
        // from-scratch recompute — i.e. the inherited word-hash and the
        // incrementally-maintained phase/loss hash never drift.
        use crate::storage::{fingerprint, phase_loss_hash, word_fingerprint};
        let mut tab = make(3);
        tab.h(0);
        tab.cnot(0, 1);
        tab.loss_channel(0, 0.3);
        tab.depolarize(1, 0.3);

        assert!(!tab.entries.dirty);
        let n = tab.entries.entries.len();
        assert_eq!(tab.entries.fingerprints.len(), n);
        assert_eq!(tab.entries.word_fingerprints.len(), n);
        assert_eq!(tab.entries.phase_loss_hashes.len(), n);
        for (i, (t, _)) in tab.entries.entries.iter().enumerate() {
            assert_eq!(
                tab.entries.word_fingerprints[i],
                word_fingerprint(t),
                "word-fingerprint cache drifted at entry {i}"
            );
            assert_eq!(
                tab.entries.phase_loss_hashes[i],
                phase_loss_hash(t),
                "phase/loss-hash cache drifted at entry {i}"
            );
            assert_eq!(
                tab.entries.fingerprints[i],
                fingerprint(t),
                "full-fingerprint cache drifted at entry {i}"
            );
        }
    }

    #[test]
    fn test_sample_shots_returns_correct_count() {
        let mut tab = make(2);
        tab.h(0);
        tab.cnot(0, 1);
        let mut sampler = tab.sampler();
        let shots = sampler.sample_shots(37);
        assert_eq!(shots.len(), 37);
        for s in shots {
            assert_eq!(s.len(), 2);
        }
    }
}
