use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::tableau::GeneralizedTableau;
use num::complex::{Complex, Complex64};
use num::traits::{One, Zero};

pub trait GeneralizedTableauTGate {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
    fn t_or_t_adj(&mut self, addr0: usize, adjoint: bool);
    fn compute_shift_z(&self, addr0: usize) -> usize;
    fn compute_phase_z(&self, addr0: usize, branch_index: usize) -> u8;
}

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

impl<const N: usize, T, C> GeneralizedTableauTGate for GeneralizedTableau<N, T, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>>,
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>> + From<Complex64>,
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

        let index_shift = self.compute_shift_z(addr0);

        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ index_shift;
            let branch_phase = self.compute_phase_z(addr0, branch_index);

            let mut phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            if adjoint {
                phase_factor.im = -phase_factor.im;
            }

            let branch_coefficient = phase_factor * coeff.clone() * complex_sin.clone();

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

    fn compute_shift_z(&self, addr0: usize) -> usize {
        // NOTE: we use LSB ordering
        let mut shift = 0usize;
        for (i, stab) in self.tableau.stabilizers.iter().enumerate() {
            shift |= (stab.word.xbits[addr0] as usize) << i;
        }
        shift
    }

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    fn compute_phase_z(&self, addr0: usize, basis_index: usize) -> u8 {
        // phase convention: 0: +1, 1: +i, 2: -1, 3: -i
        let mut phase = 0u8;
        for (i, destab) in self.tableau.destabilizers.iter().enumerate() {
            if basis_index & (1 << i) == 0 {
                // NOTE: LSB ordering; has to be consistent with shift computation
                continue;
            }

            let has_x = destab.word.xbits[addr0];
            let has_z = destab.word.zbits[addr0];

            // need to account for destabilizer phase
            phase = (phase + destab.phase) % 4;

            if has_x && has_z {
                // Y operator contributes a phase of -i
                phase = (phase + 3) % 4;
            } else if has_x {
                // X operator contributes a phase of -1
                phase = (phase + 2) % 4;
            }
        }
        phase
    }
}
