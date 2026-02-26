use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};

use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use crate::tableau::traits::TableauIndex;
use crate::traits::*;

impl<T: Config> Depolarizing<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64>,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        debug_assert!(p >= 0.0 && p <= 1.0);
        let r = rand::random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p / 3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> Depolarizing<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64>,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        debug_assert!(p >= 0.0 && p <= 1.0);
        let r = rand::random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p / 3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }
}

impl<T: Config> PauliError<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        let r = rand::random::<f64>();
        let mut cumulative = T::Coeff::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += p_.clone();
            if cumulative > r {
                match i {
                    0 => self.x(addr0),
                    1 => self.y(addr0),
                    _ => self.z(addr0),
                }
                return;
            }
        }
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> PauliError<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
        debug_assert!(p[0].clone() + p[1].clone() + p[2].clone() - 1.0 < 1e-7);
        let r = rand::random::<f64>();
        let mut cumulative = T::Coeff::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += p_.clone();
            if cumulative > r {
                match i {
                    0 => self.x(addr0),
                    1 => self.y(addr0),
                    _ => self.z(addr0),
                }
                return;
            }
        }
    }
}

fn two_qubit_pauli_error_impl<T: Config>(
    this: &mut impl Clifford,
    addr0: usize,
    addr1: usize,
    p: [T::Coeff; 15],
) where
    T::Coeff: PartialOrd<f64> + Zero,
{
    debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
    // debug_assert!(p.iter().sum() - 1.0 < 1e-7);
    let r = rand::random::<f64>();
    let sum = T::Coeff::zero();
    let idx = p
        .iter()
        .scan(sum, |acc, p_| {
            *acc += p_.clone();
            Some(acc.clone())
        })
        .position(|cum_prob| cum_prob > r);

    if let Some(i) = idx {
        #[rustfmt::skip]
        const PAULI_PAIRS: [(u8, u8); 16] = [
            (0,0),(0,1),(0,2),(0,3),
            (1,0),(1,1),(1,2),(1,3),
            (2,0),(2,1),(2,2),(2,3),
            (3,0),(3,1),(3,2),(3,3),
        ];
        let cartesian_index = PAULI_PAIRS[i + 1]; // skip II entry

        match cartesian_index.0 {
            0 => {}
            1 => this.x(addr0),
            2 => this.y(addr0),
            _ => this.z(addr0),
        }

        match cartesian_index.1 {
            0 => {}
            1 => this.x(addr1),
            2 => this.y(addr1),
            _ => this.z(addr1),
        }
    }
}

impl<T: Config> TwoQubitPauliError<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> TwoQubitPauliError<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        if self.is_lost[addr0] && self.is_lost[addr1] {
            return;
        }

        if self.is_lost[addr0] {
            // marginalize over addr0: sum columns (k fixed, j varies)
            let p_x = p[0].clone() + p[4].clone() + p[8].clone() + p[12].clone();
            let p_y = p[1].clone() + p[5].clone() + p[9].clone() + p[13].clone();
            let p_z = p[2].clone() + p[6].clone() + p[10].clone() + p[14].clone();
            self.pauli_error(addr1, [p_x, p_y, p_z]);
            return;
        }

        if self.is_lost[addr1] {
            // marginalize over addr1: sum rows (j fixed, k varies)
            let p_x = p[3].clone() + p[4].clone() + p[5].clone() + p[6].clone();
            let p_y = p[7].clone() + p[8].clone() + p[9].clone() + p[10].clone();
            let p_z = p[11].clone() + p[12].clone() + p[13].clone() + p[14].clone();
            self.pauli_error(addr0, [p_x, p_y, p_z]);
            return;
        }

        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p);
    }
}

impl<T: Config> Depolarizing2<T> for Tableau<T>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: <T as Config>::Coeff) {
        let p_arr: [T::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p_arr);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> Depolarizing2<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64> + Zero,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: <T as Config>::Coeff) {
        if self.is_lost[addr0] && self.is_lost[addr1] {
            return;
        }

        if self.is_lost[addr0] {
            self.depolarize(addr1, p.clone() * (4.0 / 5.0));
            return;
        }

        if self.is_lost[addr1] {
            self.depolarize(addr0, p * (4.0 / 5.0));
            return;
        }

        let p_arr: [T::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        two_qubit_pauli_error_impl::<T>(self, addr0, addr1, p_arr);
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> LossChannel<T>
    for GeneralizedTableau<T, I, C>
where
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64> + One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
{
    fn loss_channel(&mut self, addr0: usize, p: <T as Config>::Coeff) {
        if p < rand::random::<f64>() {
            return;
        }

        // NOTE: this is O(n^2) but also potentially removes coefficients, which is nice
        let outcome = self.measure(addr0);
        if outcome {
            // flip back to 0
            self.x(addr0);
        }
        self.is_lost[addr0] = true;
    }
}
