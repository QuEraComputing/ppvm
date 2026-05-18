use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};
use ppvm_tableau::{
    data::GeneralizedTableau, prelude::Config, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use crate::sampler::Sampler;

#[derive(Clone)]
pub struct GeneralizedTableauSum<
    T: Config,
    I,
    C: SparseVector<Complex<T::Coeff>, I> = Vec<(Complex64, I)>,
> {
    pub n_qubits: usize,
    pub(crate) entries: Vec<(GeneralizedTableau<T, I, C>, T::Coeff)>,
    pub(crate) rng: SmallRng,
    pub(crate) sum_cutoff: T::Coeff,
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableauSum<T, I, C>
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
        Self {
            n_qubits: n_qubits,
            entries: [(g_tab, T::Coeff::one())].to_vec(),
            rng: rng,
            sum_cutoff: sum_cutoff,
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
        Self {
            n_qubits: n_qubits,
            entries: [(g_tab, T::Coeff::one())].to_vec(),
            rng: rng,
            sum_cutoff: sum_cutoff,
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub(crate) fn insert_or_update_batch(
        &mut self,
        branches: &Vec<(GeneralizedTableau<T, I, C>, T::Coeff)>,
    ) {
        let mut needs_renormalize = false;
        for branch in branches.iter() {
            let dropped_any = self.insert_or_update(&branch.0, &branch.1);
            needs_renormalize |= dropped_any;
        }
        if needs_renormalize {
            self.normalize_probabilities();
        }
    }

    fn insert_or_update(&mut self, tab: &GeneralizedTableau<T, I, C>, p: &T::Coeff) -> bool {
        let idx = self
            .entries
            .iter()
            .position(|entry| Self::unsafe_equal_tableau_data_and_is_lost(&entry.0, &tab));
        let mut needs_normalize = false;
        match idx {
            Some(i) => {
                let p0 = &self.entries[i].1;
                self.entries[i].1 = p0.clone() + p.clone();
            }

            // TODO: avoid cloning here
            None => {
                if p > &self.sum_cutoff {
                    self.entries.push((tab.clone(), p.clone()));
                } else {
                    needs_normalize = true;
                }
            }
        }

        needs_normalize
    }

    pub fn truncate(&mut self) {
        let length_before_truncation = self.entries.len();
        self.entries.retain(|entry| entry.1 > self.sum_cutoff);
        if self.entries.len() < length_before_truncation {
            self.normalize_probabilities();
        }
    }

    pub fn normalize_probabilities(&mut self) {
        let norm = self
            .entries
            .iter()
            .fold(T::Coeff::zero(), |acc, entry| acc + entry.1.clone());
        for (_, p) in self.entries.iter_mut() {
            *p = p.clone() / norm.clone();
        }
    }

    fn unsafe_equal_tableau_data_and_is_lost(
        tab0: &GeneralizedTableau<T, I, C>,
        tab1: &GeneralizedTableau<T, I, C>,
    ) -> bool {
        if tab0.is_lost != tab1.is_lost {
            return false;
        }

        for (row0, row1) in tab0.tableau.data.iter().zip(tab1.tableau.data.iter()) {
            if row0.phase != row1.phase || row0.word != row1.word {
                return false;
            }
        }

        true
    }

    pub fn sampler(&mut self) -> Sampler<T, I, C> {
        self.entries
            .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut p_acc = T::Coeff::zero();
        let mut p_cum = Vec::<T::Coeff>::new();
        for entry in self.entries.iter() {
            p_acc += entry.1.clone();
            p_cum.push(p_acc.clone())
        }

        Sampler {
            p_cumulative: p_cum,
            generalized_tableau_sum: self.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::fxhash::ByteF64;
    use ppvm_runtime::traits::{Clifford, Depolarizing, LossChannel, TGate};

    type TestConfig = ByteF64<1>;
    type TestSum = GeneralizedTableauSum<TestConfig, u128>;

    fn make(n: usize) -> TestSum {
        GeneralizedTableauSum::new_with_seed(n, 1e-12, 1e-10, 42)
    }

    fn sum_of_probabilities(tab: &TestSum) -> f64 {
        tab.entries.iter().map(|e| e.1).sum()
    }

    #[test]
    fn test_new_initial_state() {
        let tab = make(3);
        assert_eq!(tab.len(), 1);
        assert_eq!(tab.n_qubits, 3);
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries[0].0.is_lost.iter().all(|x| !x));
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
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_t_gate_doesnt_add_branches() {
        let mut tab = make(2);
        tab.h(0);
        tab.t(0);
        tab.t_adj(0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_probabilities() {
        let mut tab = make(2);
        let cloned = tab.entries[0].0.clone();
        tab.entries.push((cloned, 3.0));
        tab.normalize_probabilities();
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
        assert!((tab.entries[0].1 - 0.25).abs() < 1e-12);
        assert!((tab.entries[1].1 - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_truncate_removes_low_probability_entries() {
        // sum_cutoff = 0.1
        let mut tab: TestSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.1, 42);
        let cloned = tab.entries[0].0.clone();
        tab.entries.push((cloned, 0.05));
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
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(!tab.entries[0].0.is_lost[0]);
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
            assert!((entry.1 - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn test_depolarize_zero_probability_doesnt_branch() {
        let mut tab = make(2);
        tab.depolarize(0, 0.0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
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
        assert!(tab.len() > 1, "depolarize should create branches, got len={}", tab.len());
        assert!((sum_of_probabilities(&tab) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_loss_channel_full_loss_marks_qubit_lost() {
        // p=1.0 means certain loss: original branch's probability collapses to 0 and
        // is dropped by truncate; the lost branch is renormalized to 1.0.
        let mut tab = make(2);
        tab.loss_channel(0, 1.0);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries[0].0.is_lost[0]);
        assert!(!tab.entries[0].0.is_lost[1]);
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
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
        assert!(tab.entries[0].0.is_lost[0]);
    }

    #[test]
    fn test_depolarize_skips_already_lost_qubit() {
        // Depolarizing a lost qubit must not create branches.
        let mut tab = make(2);
        tab.loss_channel(0, 1.0);
        tab.depolarize(0, 0.3);
        assert_eq!(tab.len(), 1);
        assert!((tab.entries[0].1 - 1.0).abs() < 1e-12);
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
