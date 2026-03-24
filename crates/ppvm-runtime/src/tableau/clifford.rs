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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::fxhash::ByteF64;
    use crate::tableau::CliffordExtensions;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    /// Returns (xbit, zbit, phase) for each tableau row: (destabilizer, stabilizer).
    fn rows(tab: &TestTableau) -> [(bool, bool, u8); 2] {
        [0, 1].map(|i| {
            let pw = &tab.tableau.data[i];
            (pw.word.xbits[0], pw.word.zbits[0], pw.phase)
        })
    }

    // Initial |0⟩: destabilizer = X (1,0,0), stabilizer = Z (0,1,0)

    #[test]
    fn test_sqrt_x_stabilizer() {
        // Z → -Y: forward prop √X P √X†
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x(0);
        let r = rows(&tab);
        assert_eq!(r[0], (true, false, 0), "destabilizer X should stay X");
        assert_eq!(r[1], (true, true, 2), "stabilizer Z should become -Y");
    }

    #[test]
    fn test_sqrt_x_adj_stabilizer() {
        // Z → +Y
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x_adj(0);
        let r = rows(&tab);
        assert_eq!(r[0], (true, false, 0), "destabilizer X should stay X");
        assert_eq!(r[1], (true, true, 0), "stabilizer Z should become +Y");
    }

    #[test]
    fn test_sqrt_y_stabilizer() {
        // Z → +X, X → -Z
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y(0);
        let r = rows(&tab);
        assert_eq!(r[0], (false, true, 2), "destabilizer X should become -Z");
        assert_eq!(r[1], (true, false, 0), "stabilizer Z should become +X");
    }

    #[test]
    fn test_sqrt_y_adj_stabilizer() {
        // Z → -X, X → +Z
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y_adj(0);
        let r = rows(&tab);
        assert_eq!(r[0], (false, true, 0), "destabilizer X should become +Z");
        assert_eq!(r[1], (true, false, 2), "stabilizer Z should become -X");
    }

    #[test]
    fn test_sqrt_x_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x(0);
        tab.sqrt_x_adj(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y(0);
        tab.sqrt_y_adj(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_x_fourth_power_is_identity() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        for _ in 0..4 {
            tab.sqrt_x(0);
        }
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_fourth_power_is_identity() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        for _ in 0..4 {
            tab.sqrt_y(0);
        }
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_x_on_lost_qubit_is_noop() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.is_lost[0] = true;
        tab.sqrt_x(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_on_lost_qubit_is_noop() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.is_lost[0] = true;
        tab.sqrt_y(0);
        assert_eq!(rows(&tab), initial);
    }
}
