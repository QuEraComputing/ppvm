// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Coefficient-only Heisenberg Pauli and depolarizing channels: each term is
//! rescaled in place, so no keys are created, removed, or merged. A Pauli keeps
//! its sign under errors it commutes with, so its factor is one minus twice
//! the anticommuting probability. Lost sites are skipped, matching the
//! tableau. Probabilities are assumed valid; coefficients may be complex.

use num::One;
use ppvm_traits_2::{
    Depolarizing, Depolarizing2, PauliBits, PauliError, TwoQubitPauliError, loss::LossInfo,
};

use crate::{Sum, UpdateCoefficients};

/// Packed code at `qubit` (`0=I, 1=X, 2=Z, 3=Y`), treating a lost site as `I`.
#[inline]
fn live_code<K: PauliBits + LossInfo>(key: &K, qubit: usize) -> usize {
    if key.is_lost(qubit) {
        0
    } else {
        key.pauli_code(qubit) as usize
    }
}

impl<S: UpdateCoefficients, P> PauliError<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
{
    fn pauli_error(&mut self, qubit: usize, [px, py, pz]: [S::Coeff; 3]) {
        self.check_qubit(qubit);
        let factor = |a: S::Coeff, b: S::Coeff| {
            let p = a + b;
            S::Coeff::one() - (p.clone() + p)
        };
        // Indexed by packed code: X anticommutes with Y, Z; Z with X, Y; Y with X, Z.
        let factors = [
            S::Coeff::one(),
            factor(py.clone(), pz.clone()),
            factor(px.clone(), py),
            factor(px, pz),
        ];
        self.storage.update_coefficients(|key, coeff| {
            let code = live_code(key, qubit);
            if code != 0 {
                *coeff *= &factors[code];
            }
        });
    }
}

impl<S: UpdateCoefficients, P> TwoQubitPauliError<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
{
    fn two_qubit_pauli_error(&mut self, a: usize, b: usize, p: [S::Coeff; 15]) {
        self.check_pair(a, b);
        // Probabilities use I, X, Y, Z order; packed codes use I, X, Z, Y.
        const CODE: [usize; 4] = [0, 1, 3, 2];
        let anti = |u: usize, v: usize| ((u & 1) & (v >> 1)) ^ ((u >> 1) & (v & 1)) != 0;
        // Indexed by 4 * code_a + code_b.
        let factors: [S::Coeff; 16] = std::array::from_fn(|pair| {
            let mut factor = S::Coeff::one();
            for (i, probability) in p.iter().enumerate() {
                let error = i + 1;
                if anti(pair / 4, CODE[error / 4]) != anti(pair % 4, CODE[error % 4]) {
                    factor = factor - (probability.clone() + probability.clone());
                }
            }
            factor
        });
        self.storage.update_coefficients(|key, coeff| {
            let pair = 4 * live_code(key, a) + live_code(key, b);
            if pair != 0 {
                *coeff *= &factors[pair];
            }
        });
    }
}

/// Stim's `DEPOLARIZE1(p)`: each of `X`, `Y`, `Z` with probability `p/3`.
impl<S: UpdateCoefficients, P> Depolarizing<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: From<f64>,
{
    fn depolarize1(&mut self, qubit: usize, p: S::Coeff) {
        self.check_qubit(qubit);
        let factor = S::Coeff::one() - p * S::Coeff::from(4.0 / 3.0);
        self.storage.update_coefficients(|key, coeff| {
            if live_code(key, qubit) != 0 {
                *coeff *= &factor;
            }
        });
    }
}

/// Stim's `DEPOLARIZE2(p)`: each non-identity pair with probability `p/15`.
impl<S: UpdateCoefficients, P> Depolarizing2<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: From<f64>,
{
    fn depolarize2(&mut self, a: usize, b: usize, p: S::Coeff) {
        self.check_pair(a, b);
        let factor = S::Coeff::one() - p * S::Coeff::from(16.0 / 15.0);
        self.storage.update_coefficients(|key, coeff| {
            let active_a = !key.is_lost(a) && key.pauli_code(a) != 0;
            let active_b = !key.is_lost(b) && key.pauli_code(b) != 0;
            if active_a || active_b {
                *coeff *= &factor;
            }
        });
    }
}
