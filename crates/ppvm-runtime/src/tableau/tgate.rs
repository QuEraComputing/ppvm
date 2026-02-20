use super::sparsevec::SparseVector;
use super::traits::TGate;
use crate::config::Config;
use crate::tableau::GeneralizedTableau;
use crate::traits::Coefficient;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, Zero};
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

const COS_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.8535533905932737,
    im: 0.3535533905932738,
}; // exp(im * pi / 8) * cos(pi/8)
const ISIN_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.14644660940672624,
    im: -0.3535533905932738,
}; // -im * exp(im * pi / 8) * sin(pi/8)

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },  // +1
    Complex64 { re: 0.0, im: 1.0 },  // +i
    Complex64 { re: -1.0, im: 0.0 }, // -1
    Complex64 { re: 0.0, im: -1.0 }, // -i
];

impl<T, I, C> TGate<T> for GeneralizedTableau<T, I, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I>,
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize>
        + BitOrAssign<<I as Shl<usize>>::Output>
        + BitAnd<<I as Shl<usize>>::Output, Output = I>
        + BitXor<Output = I>,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn t(&mut self, index: usize) {
        self.t_or_t_adj(index, false);
    }

    fn t_adj(&mut self, index: usize) {
        self.t_or_t_adj(index, true);
    }

    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool) {
        if self.is_lost[addr0] {
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

        self.branch_z_with_coefficients(addr0, complex_cos, complex_sin);
    }

    fn rz(&mut self, addr0: usize, theta: T::Coeff) {
        let (cos, sin) = theta.sin_cos();

        let complex_cos: Complex<T::Coeff> = Complex {
            re: cos,
            im: T::Coeff::zero(),
        };

        let i_complex_sin: Complex<T::Coeff> = Complex {
            re: T::Coeff::zero(),
            im: -sin,
        };

        self.branch_z_with_coefficients(addr0, complex_cos, i_complex_sin);
    }

    fn branch_z_with_coefficients(
        &mut self,
        addr0: usize,
        complex_cos: Complex<T::Coeff>,
        complex_sin: Complex<T::Coeff>,
    ) {
        if self.is_lost[addr0] {
            return;
        }

        let index_shift = self.compute_shift(addr0, (false, true));
        let phase_decomp = self.compute_decomposition_phase(addr0, crate::char::Pauli::Z);

        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::new();
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ index_shift;

            // get the phase contributions from duplicate destabilizers
            // and anti-commuting through destabilizers
            let branch_phase_contribution = self.compute_phase_z(addr0, idx, index_shift);
            let branch_phase = (branch_phase_contribution + phase_decomp) % 4;

            let phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            // TODO: check if we need the commented out below; isn't this accounted for by the coefficient already?
            // if adjoint {
            //     phase_factor.im = -phase_factor.im;
            // }

            let branch_coefficient = phase_factor * coeff.clone() * complex_sin.clone();
            let nonbranch_coefficient = coeff * complex_cos.clone();

            *new_coefficients
                .entry(branch_index)
                .or_insert(Complex::zero()) += branch_coefficient;
            *new_coefficients.entry(idx).or_insert(Complex::zero()) += nonbranch_coefficient;
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
