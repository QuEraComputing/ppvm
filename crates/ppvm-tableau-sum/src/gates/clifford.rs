use num::Complex;
use ppvm_runtime::{
    config::Config,
    traits::{Clifford, CliffordExtensions},
};
use ppvm_tableau::sparsevec::SparseVector;

use crate::data::GeneralizedTableauSum;
use crate::storage::EntryStore;

macro_rules! impl_generalized_tableau_sum_gate {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.entries.for_each_mut(|tab, _p| {
                tab.$name($($index), *);
            });
            // The gate mutates every entry's tableau (or no-ops on a
            // lost qubit, in which case the cached fp is still valid).
            // Conservatively clear all cached fingerprints; they'll be
            // recomputed lazily on the next insert_or_update_batch.
            self.entries.mark_keys_dirty();
        }
    };
}
pub(crate) use impl_generalized_tableau_sum_gate;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> Clifford
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_sum_gate!(x, index);
    impl_generalized_tableau_sum_gate!(y, index);
    impl_generalized_tableau_sum_gate!(z, index);
    impl_generalized_tableau_sum_gate!(h, index);
    impl_generalized_tableau_sum_gate!(s, index);
    impl_generalized_tableau_sum_gate!(cnot, control, target);
    impl_generalized_tableau_sum_gate!(cz, control, target);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> CliffordExtensions
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_sum_gate!(s_adj, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_x, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_x_adj, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_y, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_y_adj, addr0);
    impl_generalized_tableau_sum_gate!(cy, addr0, addr1);
}
