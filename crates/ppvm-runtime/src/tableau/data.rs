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

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },  // +1
    Complex64 { re: 0.0, im: 1.0 },  // +i
    Complex64 { re: -1.0, im: 0.0 }, // -1
    Complex64 { re: 0.0, im: -1.0 }, // -i
];

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
