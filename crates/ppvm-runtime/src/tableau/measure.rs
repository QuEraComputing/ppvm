use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use crate::tableau::traits::TableauIndex;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use std::ops::{BitAnd, Shl};

impl<T: Config> Measure for Tableau<T> {
    /// Measure qubit `addr0` in Z basis
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_z_anticommuting_stabilizer(addr0);
        match q {
            Some(q_idx) => {
                // Case a: random measurement outcome
                // At least one stabilizer anticommutes with Z_addr0

                // Generate random measurement outcome (50/50)
                let outcome = rand::random::<bool>();

                self.update_tableau_according_to_outcome(addr0, q_idx, outcome);

                outcome
            }
            None => {
                // Case b: deterministic measurement outcome

                self.get_deterministic_outcome(addr0)
            }
        }
    }
}

// const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
//     Complex64::new(1.0, 0.0),  // +1
//     Complex64::new(0.0, 1.0),  // +i
//     Complex64::new(-1.0, 0.0), // -1
//     Complex64::new(0.0, -1.0), // -i
// ];

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Measure for GeneralizedTableau<T, I, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    I: TableauIndex,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn measure(&mut self, addr0: usize) -> bool {
        if self.is_lost[addr0] {
            return false;
        }

        let norm: T::Coeff = self
            .coefficients
            .clone()
            .into_iter()
            .fold(T::Coeff::zero(), |acc, (v, _)| acc + v.norm_sqr());
        println!("Current norm: {:?}", norm);
        println!("{:?}", self.coefficients);

        let one_half = Complex {
            re: T::Coeff::one() * 0.5,
            im: T::Coeff::zero(),
        };
        let mut state0 = self.clone();
        state0.branch_with_coefficients(addr0, crate::char::Pauli::Z, one_half, one_half);
        let mut state1 = self.clone();
        state1.branch_with_coefficients(addr0, crate::char::Pauli::Z, one_half, -one_half);

        let prob_0 = state0
            .coefficients
            .clone()
            .into_iter()
            .fold(T::Coeff::zero(), |acc, (v, _)| acc + v.norm_sqr());
        let prob_1 = state1
            .coefficients
            .clone()
            .into_iter()
            .fold(T::Coeff::zero(), |acc, (v, _)| acc + v.norm_sqr());

        println!("p0: {:?}", prob_0.clone());
        println!("p1: {:?}", prob_1.clone());
        debug_assert!(
            (prob_0.clone() + prob_1.clone() - T::Coeff::one()) < 1e-7
                && (prob_0.clone() + prob_1.clone() - T::Coeff::one() > -1e-7)
        );

        let outcome = prob_0 < rand::random::<f64>();

        let q = self.tableau.find_z_anticommuting_stabilizer(addr0);

        match q {
            Some(q_idx) => {
                if outcome {
                    self.coefficients = state1.coefficients;
                } else {
                    self.coefficients = state0.coefficients;
                }

                self.coefficients.trim(Complex {
                    re: self.coefficient_threshold.clone(),
                    im: T::Coeff::zero(),
                });

                println!("{}", self.coefficients.len());

                self.coefficients.normalize();

                self.tableau
                    .update_tableau_according_to_outcome(addr0, q_idx, outcome);

                // let z_stabilizer_phase = self.tableau.get_deterministic_outcome(addr0);
            }
            None => {}
        };
        let z_stabilizer_phase = self.tableau.get_deterministic_outcome(addr0);

        println!("Deterministic outcome: {}", z_stabilizer_phase);
        println!("Actual outcome: {}", outcome);

        self.trim_coefficients_for_measurement(addr0, outcome ^ z_stabilizer_phase);

        outcome
    }
}
