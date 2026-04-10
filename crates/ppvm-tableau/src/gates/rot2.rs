use std::collections::HashMap;

use bitvec::view::BitView;
use num::PrimInt;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};

use crate::prelude::*;

const PAULIS: [Pauli; 4] = [Pauli::I, Pauli::X, Pauli::Z, Pauli::Y];

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> RotationTwo<T>
    for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    T::Coeff: Zero + One,
    I: TableauIndex,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
{
    fn rotate_2(
        &mut self,
        axis_a: [u8; 2],
        axis_b: [u8; 2],
        a: usize,
        b: usize,
        theta: <T as Config>::Coeff,
    ) {
        let [axis_a_x, axis_a_z] = axis_a;
        let [axis_b_x, axis_b_z] = axis_b;
        let pauli_a = PAULIS[(axis_a_z << 1 | axis_a_x) as usize];
        let pauli_b = PAULIS[(axis_b_z << 1 | axis_b_x) as usize];
        // NOTE: if both qubits are lost, the rot1 will be a no-op
        if self.is_lost[a] {
            return self.rotate_1(pauli_b, b, theta);
        } else if self.is_lost[b] {
            return self.rotate_1(pauli_a, a, theta);
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

        let mut branch_coefficients = self.coefficients.clone();
        self.compute_coefficients_after_pauli_apply(&mut branch_coefficients, b, pauli_b);
        self.compute_coefficients_after_pauli_apply(&mut branch_coefficients, a, pauli_a);

        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::new();

        for (coeff, idx) in old_coefficients {
            *new_coefficients.entry(idx).or_insert(Complex::zero()) += complex_cos * coeff;
        }

        for (branch_coeff, idx) in branch_coefficients {
            *new_coefficients.entry(idx).or_insert(Complex::zero()) += i_complex_sin * branch_coeff;
        }

        let cutoff = Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        };

        for (idx, coeff) in new_coefficients {
            if coeff.abs() > cutoff.abs() {
                self.coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::fxhash::ByteF64;
    use std::f64::consts::{FRAC_PI_2, PI};

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    // --- rxx ---

    /// rxx(π) = exp(-i·π/2·XX): Clifford, no branching, maps |00⟩ → -i|11⟩.
    #[test]
    fn test_rxx_pi_flips_both_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rxx(0, 1, PI);
        assert_eq!(tab.coefficients.len(), 1, "rxx(π) should not branch");
        assert!(tab.measure(0).unwrap());
        assert!(tab.measure(1).unwrap());
    }

    /// rxx(π/2) on |00⟩ creates a Bell-like state (|00⟩ - i|11⟩)/√2: two branches.
    #[test]
    fn test_rxx_half_pi_branches() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rxx(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "rxx(π/2) should create 2 branches"
        );
    }

    /// Two rxx(π/2) compose to rxx(π): |00⟩ → -i|11⟩ (deterministic).
    #[test]
    fn test_rxx_half_pi_twice_flips_both() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rxx(0, 1, FRAC_PI_2);
        tab.rxx(0, 1, FRAC_PI_2);
        assert!(tab.measure(0).unwrap());
        assert!(tab.measure(1).unwrap());
    }

    // --- ryy ---

    /// ryy(π) = exp(-i·π/2·YY): Clifford, no branching.
    /// YY|00⟩ = (i|1⟩)⊗(i|1⟩) = -|11⟩, so ryy(π)|00⟩ = i|11⟩.
    #[test]
    fn test_ryy_pi_flips_both_to_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.ryy(0, 1, PI);
        assert_eq!(tab.coefficients.len(), 1, "ryy(π) should not branch");
        assert!(tab.measure(0).unwrap());
        assert!(tab.measure(1).unwrap());
    }

    /// ryy(π/2) on |00⟩ creates a Bell-like state (|00⟩ + i|11⟩)/√2: two branches.
    #[test]
    fn test_ryy_half_pi_branches() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.ryy(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "ryy(π/2) should create 2 branches"
        );
    }

    /// Two ryy(π/2) compose to ryy(π): |00⟩ → i|11⟩ (deterministic).
    #[test]
    fn test_ryy_half_pi_twice_flips_both() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.ryy(0, 1, FRAC_PI_2);
        tab.ryy(0, 1, FRAC_PI_2);
        assert!(tab.measure(0).unwrap());
        assert!(tab.measure(1).unwrap());
    }

    // --- rzz ---

    /// rzz(π) on |00⟩: ZZ|00⟩ = |00⟩ (eigenvalue +1), so rzz(π)|00⟩ = -i|00⟩.
    /// No branching since ZZ is already a stabilizer of |00⟩.
    #[test]
    fn test_rzz_pi_leaves_00() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rzz(0, 1, PI);
        assert_eq!(
            tab.coefficients.len(),
            1,
            "rzz(π) on |00⟩ should not branch"
        );
        assert!(!tab.measure(0).unwrap());
        assert!(!tab.measure(1).unwrap());
    }

    /// rzz never branches on a computational basis state: ZZ is diagonal in the Z basis.
    #[test]
    fn test_rzz_does_not_branch_on_comp_basis() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rzz(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            1,
            "rzz(π/2) on |00⟩ should not branch: ZZ is a stabilizer"
        );
    }

    /// rzz(π) on |10⟩: ZZ|10⟩ = Z|1⟩⊗Z|0⟩ = -|10⟩ (eigenvalue -1), so rzz(π)|10⟩ = i|10⟩.
    #[test]
    fn test_rzz_pi_on_10() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.x(0);
        tab.rzz(0, 1, PI);
        assert_eq!(
            tab.coefficients.len(),
            1,
            "rzz(π) on |10⟩ should not branch"
        );
        assert!(tab.measure(0).unwrap());
        assert!(!tab.measure(1).unwrap());
    }

    // --- branching + physical correlation checks ---

    /// rxx(π/2)|00⟩ = (|00⟩ - i|11⟩)/√2: XX only connects |00⟩↔|11⟩, so both
    /// qubits always give the same measurement outcome.
    #[test]
    fn test_rxx_half_pi_correlated_measurements() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.rxx(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "rxx(π/2) should create 2 branches"
        );
        let r0 = tab.measure(0);
        let r1 = tab.measure(1);
        assert_eq!(
            r0, r1,
            "rxx(π/2)|00⟩ = (|00⟩-i|11⟩)/√2: measurements must agree"
        );
    }

    /// rxx(π/2)|10⟩ = (|10⟩ - i|01⟩)/√2: XX connects |10⟩↔|01⟩, so the two
    /// qubits always give opposite measurement outcomes.
    #[test]
    fn test_rxx_half_pi_on_10_anticorrelated_measurements() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.x(0);
        tab.rxx(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "rxx(π/2) on |10⟩ should create 2 branches"
        );
        let r0 = tab.measure(0);
        let r1 = tab.measure(1);
        assert_ne!(
            r0, r1,
            "rxx(π/2)|10⟩ = (|10⟩-i|01⟩)/√2: measurements must differ"
        );
    }

    /// ryy(π/2)|00⟩ = (|00⟩ + i|11⟩)/√2: YY|00⟩ = -|11⟩, so ryy(π/2) only
    /// populates |00⟩ and |11⟩. Both qubits always give the same outcome.
    #[test]
    fn test_ryy_half_pi_correlated_measurements() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.ryy(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "ryy(π/2) should create 2 branches"
        );
        let r0 = tab.measure(0);
        let r1 = tab.measure(1);
        assert_eq!(
            r0, r1,
            "ryy(π/2)|00⟩ = (|00⟩+i|11⟩)/√2: measurements must agree"
        );
    }

    /// H·rzz(π/2)·H on both qubits equals rxx(π/2) (since H·Z·H = X).
    /// Applying rzz(π/2) to |++⟩ must create 2 branches; after the second H⊗H the
    /// state is (|00⟩ - i|11⟩)/√2 and both qubits must give the same outcome.
    #[test]
    fn test_rzz_half_pi_on_plus_plus_branches_and_correlates() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.h(0);
        tab.h(1);
        tab.rzz(0, 1, FRAC_PI_2);
        assert_eq!(
            tab.coefficients.len(),
            2,
            "rzz(π/2) on |++⟩ should create 2 branches"
        );
        tab.h(0);
        tab.h(1);
        // H⊗H · rzz(π/2) · H⊗H |00⟩ = rxx(π/2)|00⟩ = (|00⟩ - i|11⟩)/√2
        let r0 = tab.measure(0);
        let r1 = tab.measure(1);
        assert_eq!(r0, r1, "H·rzz(π/2)·H = rxx(π/2): measurements must agree");
    }
}
