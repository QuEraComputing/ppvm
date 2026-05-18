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
        for branch in branches.iter() {
            self.insert_or_update(&branch.0, &branch.1);
        }
    }

    fn insert_or_update(&mut self, tab: &GeneralizedTableau<T, I, C>, p: &T::Coeff) {
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

        if needs_normalize {
            self.normalize_probabilities();
        }
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
