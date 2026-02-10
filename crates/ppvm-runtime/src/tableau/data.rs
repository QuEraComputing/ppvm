use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::phase::PhasedPauliWord;
use num::{
    One, Zero,
    complex::{Complex, Complex64},
};

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

// TODO: builder
pub struct GeneralizedTableau<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> {
    pub tableau: Tableau<N, T>,
    pub coefficients: C,
    pub is_lost: [bool; N],
    pub coefficient_threshold: T::Coeff,
}

const COS_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.8535533905932737,
    im: 0.3535533905932738,
}; // exp(im * pi / 8) * cos(pi/8)
const ISIN_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.14644660940672624,
    im: -0.3535533905932738,
}; // -im * exp(im * pi / 8) * sin(pi/8)

impl<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> GeneralizedTableau<N, T, C>
where
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>> + From<Complex64>,
{
    pub fn new(coefficient_threshold: T::Coeff) -> Self {
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
            coefficient_threshold,
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

        let complex_cos: Complex<T::Coeff> = if adjoint {
            COS_PI_OVER_8_TIMES_EXPIPI8.conj().into()
        } else {
            COS_PI_OVER_8_TIMES_EXPIPI8.into()
        };

        let complex_sin: Complex<T::Coeff> = if adjoint {
            ISIN_PI_OVER_8_TIMES_EXPIPI8.conj().into()
        } else {
            ISIN_PI_OVER_8_TIMES_EXPIPI8.into()
        };

        let index_shift = self.compute_shift_z(index);

        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_coefficient = coeff.clone() * complex_sin.clone();

            // TODO: phase
            let branch_index = idx ^ index_shift;

            let nonbranch_coefficient = coeff * complex_cos.clone();
            self.coefficients
                .add_or_insert(branch_index, branch_coefficient);
            self.coefficients.add_or_insert(idx, nonbranch_coefficient);
        }

        // TODO: more efficient trimming above
        self.coefficients.trim(Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        });
    }

    fn compute_shift_z(&self, index: usize) -> usize {
        let mut shift = 0usize;
        for (i, stab) in self.tableau.stabilizers.iter().enumerate() {
            shift |= (stab.word.xbits[index] as usize) << i;
        }
        shift
    }
}
