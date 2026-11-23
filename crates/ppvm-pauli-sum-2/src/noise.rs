// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Coefficient-only Heisenberg Pauli and depolarizing channels: each term is
//! rescaled in place, so no keys are created, removed, or merged. Lost sites
//! are skipped, matching the tableau, where gates and errors skip lost qubits.
//! Probabilities are assumed valid; coefficients may be complex.

use num::One;
use ppvm_traits_2::{
    Depolarizing, Depolarizing2, PauliBits, PauliError, TwoQubitPauliError, loss::LossInfo,
};

use crate::{Sum, UpdateCoefficients};

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
        let factors = [
            S::Coeff::one(),
            factor(py.clone(), pz.clone()),
            factor(px.clone(), py),
            factor(px, pz),
        ];
        self.storage.update_coefficients(|key, coeff| {
            if !key.is_lost(qubit) && key.pauli_code(qubit) != 0 {
                *coeff *= &factors[key.pauli_code(qubit) as usize];
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
        // Probabilities use I,X,Y,Z order; packed word codes use I,X,Z,Y.
        let codes = [0u8, 1, 3, 2];
        let factors: [S::Coeff; 16] = std::array::from_fn(|pair| {
            let mut factor = S::Coeff::one();
            for (i, probability) in p.iter().enumerate() {
                let error = i + 1;
                let anti = |x: u8, y: u8| ((x & 1) & (y >> 1)) ^ ((x >> 1) & (y & 1));
                if anti((pair / 4) as u8, codes[error / 4])
                    ^ anti((pair % 4) as u8, codes[error % 4])
                    != 0
                {
                    factor = factor - (probability.clone() + probability.clone());
                }
            }
            factor
        });
        self.storage.update_coefficients(|key, coeff| {
            let x = if key.is_lost(a) {
                0
            } else {
                key.pauli_code(a) as usize
            };
            let y = if key.is_lost(b) {
                0
            } else {
                key.pauli_code(b) as usize
            };
            if x != 0 || y != 0 {
                *coeff *= &factors[x * 4 + y];
            }
        });
    }
}

impl<S: UpdateCoefficients, P> Depolarizing<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: From<f64>,
{
    fn depolarize1(&mut self, qubit: usize, p: S::Coeff) {
        self.check_qubit(qubit);
        let factor = S::Coeff::one() - p * S::Coeff::from(4.0 / 3.0);
        self.storage.update_coefficients(|key, coeff| {
            if !key.is_lost(qubit) && key.pauli_code(qubit) != 0 {
                *coeff *= &factor;
            }
        });
    }
}

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
