use super::data::Tableau;
use crate::config::Config;
use crate::traits::Clifford;

macro_rules! impl_tableau_clifford {
    ($name:ident, $($index:ident),*) => {
        fn $name(&mut self, $($index: usize),*) {
            // self.destabilizers.iter_mut().for_each(|pw| {
            //     pw.$name($($index),*);
            // });
            self.stabilizers.iter_mut().for_each(|pw| {
                pw.$name($($index),*);
            });
        }
    };
}

impl<const N: usize, T: Config> Clifford for Tableau<N, T> {
    impl_tableau_clifford!(x, index);
    impl_tableau_clifford!(y, index);
    impl_tableau_clifford!(z, index);
    impl_tableau_clifford!(h, index);
    impl_tableau_clifford!(s, index);
    impl_tableau_clifford!(cnot, control, target);
    impl_tableau_clifford!(cz, control, target);
}
