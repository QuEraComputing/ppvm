use crate::char::Pauli;
use crate::loss::LossyPauliWord;
use crate::traits::*;
use crate::{config::Config, sum::PauliSum};
use std::hash::BuildHasher;

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
                Pauli::L => {}
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
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::I, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[8].clone();
            }

            (Pauli::I, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::X, Pauli::I) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::X, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone();
            }

            (Pauli::X, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::X, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[7].clone();
            }

            (Pauli::Y, Pauli::I) => {
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

            (Pauli::Y, Pauli::X) => {
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

            (Pauli::Y, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[8].clone();
            }

            (Pauli::Y, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[11].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::Z, Pauli::I) => {
                *v *= 1.0f64
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::Z, Pauli::X) => {
                *v *= 1.0f64
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[4].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[8].clone();
            }

            (Pauli::Z, Pauli::Y) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[14].clone()
                    - 2.0f64 * p[2].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[5].clone()
                    - 2.0f64 * p[7].clone()
                    - 2.0f64 * p[9].clone();
            }

            (Pauli::Z, Pauli::Z) => {
                *v *= 1.0f64
                    - 2.0f64 * p[0].clone()
                    - 2.0f64 * p[1].clone()
                    - 2.0f64 * p[10].clone()
                    - 2.0f64 * p[12].clone()
                    - 2.0f64 * p[13].clone()
                    - 2.0f64 * p[3].clone()
                    - 2.0f64 * p[6].clone()
                    - 2.0f64 * p[7].clone();
            }

            _ => {
                // TODO: no action on loss needs test!
            }
        })
    }
}

impl<T: Config> DepolarizingError<T> for PauliSum<T>
where
    f64: std::ops::Mul<T::Coeff, Output = T::Coeff>
        + std::ops::Add<T::Coeff, Output = T::Coeff>
        + std::ops::Sub<T::Coeff, Output = T::Coeff>,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        self.scale(|k, v| match k.get(addr0) {
            Pauli::I => {}
            Pauli::X => {
                *v *= 1.0f64 - 4.0f64 / 3.0f64 * p.clone();
            }
            Pauli::Y => {
                *v *= 1.0f64 - 4.0f64 / 3.0f64 * p.clone();
            }
            Pauli::Z => {
                *v *= 1.0f64 - 4.0f64 / 3.0f64 * p.clone();
            }
            Pauli::L => {}
        });
    }
}

/// Loss channel implementation for PauliSum
///
/// This trait reduces the trace of the density matrix as (1 - p) per lost qubit.
/// While this is technically correct, you may want to count loss as a contribution
/// to the zero state of a qubit. Refer to `LossyPauliWord` and the `ResetLossChannel`
/// trait for that functionality.
impl<T: Config> LossChannel<T> for PauliSum<T>
where
    f64: std::ops::Sub<T::Coeff, Output = T::Coeff>,
{
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff) {
        self.map_insert(|k, v| match k.get(addr0) {
            Pauli::L => {
                let new_v = v.clone() * p.clone();
                let mut new_k = k.clone();
                new_k.set(addr0, Pauli::I);
                Some((new_k, new_v))
            }
            Pauli::I | Pauli::X | Pauli::Y | Pauli::Z => {
                *v *= 1.0f64 - p.clone();
                None
            }
        });
    }
}

/// Reset-loss channel implementation for PauliSum.
///
/// This trait is **only implemented for `LossyPauliWord`** and cannot be used with
/// regular `PauliWord`, as the loss state cannot represent loss.
impl<S, H, T> ResetLossChannel<T> for PauliSum<T>
where
    S: PauliStorage,
    H: BuildHasher + Clone + Default,
    T: Config<PauliWordType = LossyPauliWord<S, H>>,
{
    /// Apply the reset-loss channel to qubit at `addr0`.
    fn reset_loss_channel(&mut self, addr0: usize) {
        self.map_insert(|k, v| match k.get(addr0) {
            Pauli::L => {
                *v *= 0.0;
                None
            }
            Pauli::I | Pauli::Z => {
                let mut new_k = k.clone();
                new_k.set(addr0, Pauli::L);
                Some((new_k, v.clone()))
            }
            Pauli::X | Pauli::Y => None,
        });
    }
}
