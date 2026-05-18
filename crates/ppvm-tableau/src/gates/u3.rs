// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use crate::prelude::*;
use num::Complex;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> U3Gate<T> for GeneralizedTableau<T, I, C>
where
    GeneralizedTableau<T, I, C>: RotationOne<T>,
{
    fn u3(&mut self, addr0: usize, theta: T::Coeff, phi: T::Coeff, lambda: T::Coeff) {
        // RZ(phi)RY(theta)RZ(lambda)
        self.rz(addr0, lambda);
        self.ry(addr0, theta);
        self.rz(addr0, phi);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::fxhash::ByteF64;
    use std::f64::consts::PI;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    #[test]
    fn test_u3_identity() {
        // U3(0, 0, 0) = I, |0⟩ stays |0⟩
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.u3(0, 0.0, 0.0, 0.0);
        assert_eq!(tab.coefficients.len(), 1);
        assert!(!tab.measure(0).unwrap());
    }

    #[test]
    fn test_u3_x_gate() {
        // U3(π, 0, 0) = RY(π) (up to global phase): |0⟩ → |1⟩
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.u3(0, PI, 0.0, 0.0);
        assert_eq!(tab.coefficients.len(), 1);
        assert!(tab.measure(0).unwrap());
    }

    #[test]
    fn test_u3_superposition() {
        // U3(π/2, 0, 0) = RY(π/2): creates a superposition
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.u3(0, PI / 2.0, 0.0, 0.0);
        assert_eq!(tab.coefficients.len(), 2);
    }

    #[test]
    fn test_u3_matches_rz_ry_rz() {
        // U3(θ, φ, λ) must give identical measurement statistics to manual RZ(φ)·RY(θ)·RZ(λ)
        let (theta, phi, lambda) = (0.34 * PI, 0.21 * PI, 0.46 * PI);

        let mut tab_u3: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab_u3.u3(0, theta, phi, lambda);

        let mut tab_manual: TestTableau = GeneralizedTableau::new_with_seed(1, 1e-12, 0);
        tab_manual.rz(0, lambda);
        tab_manual.ry(0, theta);
        tab_manual.rz(0, phi);

        for seed in 0..200 {
            let result_u3 = tab_u3.fork(Some(seed)).measure(0).unwrap();
            let result_manual = tab_manual.fork(Some(seed)).measure(0).unwrap();
            assert_eq!(result_u3, result_manual, "mismatch at seed {}", seed);
        }
    }
}
