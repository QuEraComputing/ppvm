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
    use std::f64::consts::PI;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    /// rx(π) = exp(-i·π/2·X) acts as -iX on |0⟩, giving -i|1⟩.
    #[test]
    fn test_rx_pi_flips_zero_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        println!("{}", tab);
        tab.rx(0, PI);
        println!("{}", tab);
        let result = tab.measure(0);
        assert!(result, "Expected to measure 1, got {}", result);
    }
}
