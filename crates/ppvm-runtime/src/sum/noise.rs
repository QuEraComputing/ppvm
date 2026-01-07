use crate::char::Pauli;
use crate::traits::*;
use crate::{config::Config, sum::PauliSum};

impl<T: Config> PauliError<T> for PauliSum<T>
where
    f64: std::ops::Mul<T::Coeff, Output = T::Coeff>
        + std::ops::Add<T::Coeff, Output = T::Coeff>
        + std::ops::Sub<T::Coeff, Output = T::Coeff>,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        self.scale(|k, v| {
            match k.get(addr0) {
                Pauli::I => {}
                Pauli::X => {
                    *v *= 1.0f64 - 2.0f64 * p[1].clone() - 2.0f64 * p[2].clone();
                }
                Pauli::Y => {
                    *v *= 1.0f64 - 2.0f64 * p[0].clone() - 2.0f64 * p[2].clone();
                }
                Pauli::Z => {
                    *v *= 1.0f64 - 2.0f64 * p[0].clone() - 2.0f64 * p[1].clone();
                }
            };
        });
    }
}

impl<T: Config> TwoQubitPauliError<T> for PauliSum<T>
where
    f64: std::ops::Mul<T::Coeff, Output = T::Coeff>
        + std::ops::Add<T::Coeff, Output = T::Coeff>
        + std::ops::Sub<T::Coeff, Output = T::Coeff>,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        self.scale(|k, v| match (k.get(addr0), k.get(addr1)) {
            (Pauli::I, Pauli::I) => {}
            (Pauli::I, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::I, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::I, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[8].clone();
            }
            (Pauli::X, Pauli::I) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone();
            }
            (Pauli::X, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[8].clone();
            }
            (Pauli::X, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[7].clone();
            }
            (Pauli::X, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Y, Pauli::I) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[15].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Y, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone();
            }
            (Pauli::Y, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone();
            }
            (Pauli::Y, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Z, Pauli::I) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Z, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Z, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[9].clone();
            }
            (Pauli::Z, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone();
            }
        })
    }
}
