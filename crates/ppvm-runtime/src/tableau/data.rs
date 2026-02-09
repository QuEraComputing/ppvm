use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::phase::PhasedPauliWord;
use itertools::Itertools;
use num::{Num, One, Zero, complex::Complex};

#[derive(Clone, Debug)]
pub struct Tableau<const N: usize, T: Config> {
    pub destabilizers: [PhasedPauliWord<T::Storage, T::BuildHasher>; N],
    pub stabilizers: [PhasedPauliWord<T::Storage, T::BuildHasher>; N],
}

impl<const N: usize, T: Config> Tableau<N, T> {
    pub fn new() -> Self {
        let stabilizers = std::array::from_fn(|i| {
            let mut pw = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(N);
            pw.set(i, crate::char::Pauli::Z);
            pw
        });
        let destabilizers = std::array::from_fn(|i| {
            let mut pw = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(N);
            pw.set(i, crate::char::Pauli::X);
            pw
        });
        Self {
            destabilizers,
            stabilizers,
        }
    }
}

pub struct GeneralizedTableau<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> {
    pub tableau: Tableau<N, T>,
    pub coefficients: C,
    pub is_lost: [bool; N],
}

const COS_PI_OVER_8: f64 = 0.9238795325112867; // cos(pi/8)
const SIN_PI_OVER_8: f64 = 0.3826834323650898; // sin(pi/8)

impl<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> GeneralizedTableau<N, T, C>
where
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>,
{
    pub fn new() -> Self {
        let mut coefficients = C::new();
        let complex_one = Complex {
            re: T::Coeff::one(),
            im: T::Coeff::zero(),
        };
        coefficients.unsafe_insert(0, complex_one);
        Self {
            tableau: Tableau::new(),
            coefficients: coefficients,
            is_lost: [false; N],
        }
    }

    pub fn t(&mut self, index: usize) {
        self.t_or_t_adj(index, false);
    }

    pub fn t_adj(&mut self, index: usize) {
        self.t_or_t_adj(index, true);
    }

    fn t_or_t_adj(&mut self, index: usize, adjoint: bool) {
        if self.is_lost[index] {
            return;
        }

        let complex_cos = Complex {
            re: COS_PI_OVER_8.into(),
            im: T::Coeff::zero(),
        };

        let complex_sin = if adjoint {
            Complex {
                re: T::Coeff::zero(),
                im: SIN_PI_OVER_8.into(),
            }
        } else {
            Complex {
                re: T::Coeff::zero(),
                im: (-SIN_PI_OVER_8).into(),
            }
        };

        let mut new_coefficients = C::new();
        for (coeff, idx) in self.coefficients.clone().into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_coefficient = coeff.clone() * complex_sin.clone();

            // TODO: phase
            let branch_index = self.compute_shift_z(idx);

            let nonbranch_coefficient = coeff.clone() * complex_cos.clone();
            new_coefficients.add_or_insert(branch_index, branch_coefficient);
            new_coefficients.add_or_insert(idx, nonbranch_coefficient);
        }

        self.coefficients = new_coefficients;
    }

    fn compute_shift_z(&self, index: usize) -> usize {
        let index_bits = (0..N).map(|i| ((index >> i) & 1) != 0); // LSB-first
        let beta = self.tableau.stabilizers.iter().map(|k| k.word.xbits[index]);

        // mod 2 addition
        let shift_bits = index_bits.zip(beta).map(|(i, b)| i ^ b);

        let shift = shift_bits
            .into_iter()
            .enumerate()
            .fold(0, |acc, (i, b)| acc | ((b as usize) << i));

        shift
    }
}
