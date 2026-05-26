pub mod entry_store;
pub mod vec;

pub use entry_store::EntryStore;
use fxhash::{FxHashMap, FxHasher};
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use std::{
    hash::{Hash, Hasher},
    ops::AddAssign,
};

pub(crate) fn fingerprint<T, I, C>(tab: &GeneralizedTableau<T, I, C>) -> u64
where
    T: Config,
    I:,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    let mut hasher = FxHasher::default();
    tab.is_lost.hash(&mut hasher);
    for row in tab.tableau.data.iter() {
        row.phase.hash(&mut hasher);
        row.word.hash(&mut hasher);
    }
    hasher.finish()
}

pub(crate) fn structurally_equal<T, I, C>(
    tab0: &GeneralizedTableau<T, I, C>,
    tab1: &GeneralizedTableau<T, I, C>,
) -> bool
where
    T: Config,
    T::Coeff: One + Zero + Clone + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    // NOTE: comparing is_lost and rows is only necessary to avoid hash collisions

    if tab0.is_lost != tab1.is_lost {
        return false;
    }

    if tab0.coefficients.len() != tab1.coefficients.len() {
        return false;
    }

    // Cheaper row comparison first; coefficient compare is O(K) below.
    for (row0, row1) in tab0.tableau.data.iter().zip(tab1.tableau.data.iter()) {
        if row0.phase != row1.phase || row0.word != row1.word {
            return false;
        }
    }

    // Build a HashMap over tab1's coefficients so each tab0 lookup is O(1).
    let mut tab1_map: FxHashMap<I, Complex<T::Coeff>> =
        FxHashMap::with_capacity_and_hasher(tab1.coefficients.len(), Default::default());
    for (val, idx) in tab1.coefficients.iter() {
        tab1_map.insert(*idx, *val);
    }

    let threshold_sq = tab0.coefficient_threshold.clone() * tab0.coefficient_threshold.clone();
    let zero = Complex {
        re: T::Coeff::zero(),
        im: T::Coeff::zero(),
    };
    for (val0, idx0) in tab0.coefficients.iter() {
        let val1 = tab1_map.get(idx0).copied().unwrap_or(zero);
        if (*val0 - val1).norm_sqr() >= threshold_sq {
            return false;
        }
    }

    true
}
