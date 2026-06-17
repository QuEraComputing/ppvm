// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::char::Pauli;
use crate::loss::LossyPauliWord;
use crate::traits::*;
use crate::{config::Config, sum::PauliSum};
use num::traits::Float;
use std::hash::BuildHasher;

#[inline(always)]
fn pauli_code<W: PauliWordTrait>(word: &W, addr: usize) -> usize {
    (word.get_xbit(addr) as usize) | ((word.get_zbit(addr) as usize) << 1)
}

impl<T: Config> PauliError<T> for PauliSum<T> {
    fn pauli_error(&mut self, targets: impl Targets, p: [<T as Config>::Coeff; 3]) {
        for addr0 in targets.each() {
            let p = p.clone();
            let one = T::Coeff::from(1.0);
            let x_factor = one.clone() - p[1].clone() * 2.0 - p[2].clone() * 2.0;
            let z_factor = one.clone() - p[0].clone() * 2.0 - p[1].clone() * 2.0;
            let y_factor = one - p[0].clone() * 2.0 - p[2].clone() * 2.0;

            self.scale(move |k, v| {
                if k.get_lbit(addr0) {
                    return;
                }
                match pauli_code(k, addr0) {
                    0 => {}
                    1 => *v *= x_factor.clone(),
                    2 => *v *= z_factor.clone(),
                    3 => *v *= y_factor.clone(),
                    _ => unreachable!(),
                }
            });
        }
    }
}

/// Helper: compute `1 - 2 * (p[i_0] + p[i_1] + ... + p[i_n])`.
#[inline]
fn one_minus_two_sum<C: Coefficient, const N: usize>(p: &[C; 15], indices: [usize; N]) -> C {
    let mut acc = C::from(1.0);
    for i in indices {
        acc = acc - p[i].clone() * 2.0;
    }
    acc
}

impl<T: Config> TwoQubitPauliError<T> for PauliSum<T> {
    fn two_qubit_pauli_error(&mut self, targets: impl Targets, p: [<T as Config>::Coeff; 15]) {
        for (addr0, addr1) in targets.pairs() {
        let p = p.clone();
        self.scale(|k, v| match (k.get(addr0), k.get(addr1)) {
            (Pauli::I, Pauli::I) => {}
            (Pauli::I, Pauli::X) => {
                *v *= one_minus_two_sum(&p, [1, 10, 13, 14, 2, 5, 6, 9]);
            }
            (Pauli::I, Pauli::Y) => {
                *v *= one_minus_two_sum(&p, [0, 10, 12, 14, 2, 4, 6, 8]);
            }
            (Pauli::I, Pauli::Z) => {
                *v *= one_minus_two_sum(&p, [0, 1, 12, 13, 4, 5, 8, 9]);
            }
            (Pauli::X, Pauli::I) => {
                *v *= one_minus_two_sum(&p, [10, 11, 12, 13, 14, 7, 8, 9]);
            }
            (Pauli::X, Pauli::X) => {
                *v *= one_minus_two_sum(&p, [1, 11, 12, 2, 5, 6, 7, 8]);
            }
            (Pauli::X, Pauli::Y) => {
                *v *= one_minus_two_sum(&p, [0, 11, 13, 2, 4, 6, 7, 9]);
            }
            (Pauli::X, Pauli::Z) => {
                *v *= one_minus_two_sum(&p, [0, 1, 10, 11, 14, 4, 5, 7]);
            }
            (Pauli::Y, Pauli::I) => {
                *v *= one_minus_two_sum(&p, [11, 12, 13, 14, 3, 4, 5, 6]);
            }
            (Pauli::Y, Pauli::X) => {
                *v *= one_minus_two_sum(&p, [1, 10, 11, 12, 2, 3, 4, 9]);
            }
            (Pauli::Y, Pauli::Y) => {
                *v *= one_minus_two_sum(&p, [0, 10, 11, 13, 2, 3, 5, 8]);
            }
            (Pauli::Y, Pauli::Z) => {
                *v *= one_minus_two_sum(&p, [0, 1, 11, 14, 3, 6, 8, 9]);
            }
            (Pauli::Z, Pauli::I) => {
                *v *= one_minus_two_sum(&p, [10, 3, 4, 5, 6, 7, 8, 9]);
            }
            (Pauli::Z, Pauli::X) => {
                *v *= one_minus_two_sum(&p, [1, 13, 14, 2, 3, 4, 7, 8]);
            }
            (Pauli::Z, Pauli::Y) => {
                *v *= one_minus_two_sum(&p, [0, 12, 14, 2, 3, 5, 7, 9]);
            }
            (Pauli::Z, Pauli::Z) => {
                *v *= one_minus_two_sum(&p, [0, 1, 10, 12, 13, 3, 6, 7]);
            }
            _ => {
                // NOTE: if just one atom is lost, then there is no
                // well-defined noise channel on the other atom
                // so we don't apply any noise
            }
        })
        }
    }
}

impl<T: Config> Depolarizing<T> for PauliSum<T> {
    fn depolarize1(&mut self, targets: impl Targets, p: T::Coeff) {
        for addr0 in targets.each() {
            let factor = T::Coeff::from(1.0) - p.clone() * (4.0 / 3.0);
            self.scale(move |k, v| {
                if !k.get_lbit(addr0) && pauli_code(k, addr0) != 0 {
                    *v *= factor.clone();
                }
            });
        }
    }
}

impl<T: Config> Depolarizing2<T> for PauliSum<T> {
    fn depolarize2(&mut self, targets: impl Targets, p: T::Coeff) {
        for (addr0, addr1) in targets.pairs() {
            let factor = T::Coeff::from(1.0) - p.clone() * (16.0 / 15.0);
            self.scale(move |k, v| {
                if k.get_lbit(addr0) || k.get_lbit(addr1) {
                    return;
                }
                if pauli_code(k, addr0) != 0 || pauli_code(k, addr1) != 0 {
                    *v *= factor.clone();
                }
            });
        }
    }
}

impl<T: Config> AmplitudeDamping<T> for PauliSum<T>
where
    T::Coeff: Float,
{
    fn amplitude_damping(&mut self, addr0: usize, gamma: <T as Config>::Coeff) {
        self.map_insert(|k, v| match k.get(addr0) {
            Pauli::I | Pauli::L => None,

            Pauli::X | Pauli::Y => {
                *v *= (T::Coeff::from(1.0) - gamma).sqrt();
                None
            }

            Pauli::Z => {
                // branch to gamma * I
                let new_v = *v * gamma;
                let mut new_k = k.clone();
                new_k.set(addr0, Pauli::I);

                *v *= T::Coeff::from(1.0) - gamma;

                Some((new_k, new_v))
            }
        });
    }
}

/// Loss channel implementation for PauliSum
///
/// This trait reduces the trace of the density matrix as (1 - p) per lost qubit.
/// While this is technically correct, you may want to count loss as a contribution
/// to the zero state of a qubit. Refer to `LossyPauliWord` and the `ResetLossChannel`
/// trait for that functionality.
impl<T: Config> LossChannel<T> for PauliSum<T> {
    fn loss_channel(&mut self, addr0: usize, p: T::Coeff) {
        self.map_insert(|k, v| match k.get(addr0) {
            Pauli::L => {
                let new_v = v.clone() * p.clone();
                let mut new_k = k.clone();
                new_k.set(addr0, Pauli::I);
                Some((new_k, new_v))
            }
            Pauli::I | Pauli::X | Pauli::Y | Pauli::Z => {
                *v *= T::Coeff::from(1.0) - p.clone();
                None
            }
        });
    }
}

impl<T: Config> CorrelatedLossChannel<T> for PauliSum<T> {
    /// Apply a correlated loss channel to qubits at `addr0` and `addr1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: The probability of losing both qubits simultaneously when
    ///   both of them are in the qubit subspace.
    /// * `p[1]`: The probability of losing either one qubit when both of them are
    ///   in the qubit subspace.
    /// * `p[2]`: The probability of losing one qubit when the other one has already
    ///   been lost prior to the channel.
    fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 3]) {
        self.map_insert_multiple(|k, v| {
            match (k.get(addr0), k.get(addr1)) {
                (Pauli::L, Pauli::L) => {
                    // both qubits lost
                    let v_il = v.clone() * p[2].clone();
                    let mut k_il = k.clone();
                    k_il.set(addr0, Pauli::I);
                    k_il.set(addr1, Pauli::L);
                    let mut k_li = k.clone();
                    k_li.set(addr0, Pauli::L);
                    k_li.set(addr1, Pauli::I);

                    let v_ii = v.clone() * p[0].clone();
                    let mut k_ii = k.clone();
                    k_ii.set(addr0, Pauli::I);
                    k_ii.set(addr1, Pauli::I);

                    Some(Vec::from([
                        (k_il, v_il.clone()),
                        (k_li, v_il),
                        (k_ii, v_ii),
                    ]))
                }

                (_, Pauli::L) => {
                    // case qubit 0 in qubit subspace, qubit 1 is lost
                    let mut new_k = k.clone();
                    new_k.set(addr1, Pauli::I);
                    let new_v = v.clone() * p[1].clone();

                    *v *= T::Coeff::from(1.0) - p[2].clone();

                    Some(Vec::from([(new_k, new_v)]))
                }

                (Pauli::L, _) => {
                    // case qubit 0 is lost, qubit 1 in qubit subspace

                    let mut new_k = k.clone();
                    new_k.set(addr0, Pauli::I);
                    let new_v = v.clone() * p[1].clone();

                    *v *= T::Coeff::from(1.0) - p[2].clone();

                    Some(Vec::from([(new_k, new_v)]))
                }

                (_, _) => {
                    // case both qubits in qubit subspace
                    *v *= T::Coeff::from(1.0) - p[1].clone() * 2.0 - p[0].clone();
                    None
                }
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
