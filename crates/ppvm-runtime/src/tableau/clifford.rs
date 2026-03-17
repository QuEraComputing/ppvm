use super::data::{GeneralizedTableau, Tableau};
use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::tableau::CliffordExtensions;
use crate::traits::Clifford;
use num::complex::Complex;

macro_rules! impl_tableau_clifford {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            self.data.iter_mut().for_each(|pw| {
                pw.$name($($index),*);
            });
        }
    };
}

macro_rules! impl_generalized_tableau_clifford {
    ($name:ident, $index:ident) => {
        fn $name(&mut self, $index: usize) {
            if self.is_lost[$index] {
                return;
            }
            self.tableau.$name($index);
        }
    };
    ($name:ident, $index0:ident, $index1:ident) => {
        fn $name(&mut self, $index0: usize, $index1: usize) {
            if self.is_lost[$index0] || self.is_lost[$index1] {
                return;
            }
            self.tableau.$name($index0, $index1);
        }
    };
}

impl<T: Config> Clifford for Tableau<T> {
    impl_tableau_clifford!(x, index);
    impl_tableau_clifford!(y, index);
    impl_tableau_clifford!(z, index);
    impl_tableau_clifford!(h, index);
    impl_tableau_clifford!(cnot, control, target);
    impl_tableau_clifford!(cz, control, target);

    fn s(&mut self, index: usize) {
        // NOTE: S is the only clifford where forward and backward propagation differ
        // since it's non-hermitian
        // only difference is the phase though
        // TODO: just use the conjugate sdagger impl
        self.data.iter_mut().for_each(|pw| {
            let phase = (pw.word.xbits[index] & pw.word.zbits[index]) as u8;
            pw.word.s(index);
            pw.add_phase(phase << 1);
        });
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Clifford for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_clifford!(x, index);
    impl_generalized_tableau_clifford!(y, index);
    impl_generalized_tableau_clifford!(z, index);
    impl_generalized_tableau_clifford!(h, index);
    impl_generalized_tableau_clifford!(s, index);
    impl_generalized_tableau_clifford!(cnot, control, target);
    impl_generalized_tableau_clifford!(cz, control, target);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordExtensions
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    fn s_adj(&mut self, addr0: usize) {
        if self.is_lost[addr0] {
            return;
        }

        // NOTE: the backwards prop version of S is just S_adj
        self.tableau.data.iter_mut().for_each(|pw| {
            pw.s(addr0);
        });
    }
}
