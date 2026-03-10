use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};

use crate::config::Config;
use crate::tableau::GeneralizedTableau;
use crate::tableau::sparsevec::SparseVector;
use crate::traits::*;
use crate::{char::Pauli, tableau::traits::TableauIndex};

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> RotationOne<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: Zero + One,
    I: TableauIndex,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
{
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: <T as Config>::Coeff) {
        if self.is_lost[addr0] {
            return;
        }
        let (sin, cos) = (theta * 0.5.into()).sin_cos();

        let complex_cos: Complex<T::Coeff> = Complex {
            re: cos,
            im: T::Coeff::zero(),
        };

        let i_complex_sin: Complex<T::Coeff> = Complex {
            re: T::Coeff::zero(),
            im: -sin,
        };

        self.branch_with_coefficients(addr0, axis, complex_cos, i_complex_sin);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::fxhash::ByteF64;
    use crate::tableau::Measure;
    use std::f64::consts::{FRAC_PI_2, PI};

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    /// rx(π) = exp(-i·π/2·X) acts as -iX on |0⟩, giving -i|1⟩.
    #[test]
    fn test_rx_pi_flips_zero_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.rx(0, PI);

        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);

        let result = tab.measure(0);
        assert!(result, "Expected to measure 1, got {}", result);
    }

    #[test]
    fn test_ry_pi_flips_zero_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.ry(0, PI);
        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);
        let result = tab.measure(0);
        assert!(result, "Expected to measure 1, got {}", result);
    }

    #[test]
    fn test_rz_leaves_zero_invariant() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.rz(0, 0.123);
        // make sure we don't branch
        assert_eq!(tab.coefficients.len(), 1);
        let result = tab.measure(0);
        assert!(!result, "Expected to measure 0, got {}", result);
    }

    #[test]
    fn test_rx_pi_2_statistics() {
        let mut tab: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab.rx(0, FRAC_PI_2);
        // make sure we branch
        assert_eq!(tab.coefficients.len(), 2);

        let trials = 100;
        let mut count_one = 0;
        for i in 0..trials {
            let mut tmp_tab = tab.fork(Some(i));
            let result = tmp_tab.measure(0);
            count_one += result as u8;
        }
        println!("{}", count_one);
        assert!(35 < count_one && count_one < 65);
    }

    #[test]
    fn test_ry_round_trip() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-10);
        tab.ry(0, 2.0 * PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(0));
    }

    #[test]
    fn test_two_qubit_case() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-10);
        tab.rx(0, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(1));
        assert!(tab.measure(0));

        tab.rx(1, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(tab.measure(1));
        assert!(tab.measure(0));

        tab.rx(0, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(tab.measure(1));
        assert!(!tab.measure(0));

        tab.rx(1, PI);

        assert_eq!(tab.coefficients.len(), 1);

        assert!(!tab.measure(1));
        assert!(!tab.measure(0));
    }
}
