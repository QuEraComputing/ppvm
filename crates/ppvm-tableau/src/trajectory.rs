// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Explicit stochastic-boundary helpers for trajectory caching.
//!
//! Normal noise-channel methods intentionally hide their RNG draw. Cached
//! execution needs that draw as a cache-edge key, so these helpers sample a
//! compact branch code and apply only the deterministic operation for the
//! sampled branch.

use std::fmt::Debug;

use bitvec::view::BitView;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::{One, ToPrimitive, Zero};
use rand::RngExt;

use crate::prelude::*;

const PAULI_PAIRS: [(u8, u8); 16] = [
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 0),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 0),
    (2, 1),
    (2, 2),
    (2, 3),
    (3, 0),
    (3, 1),
    (3, 2),
    (3, 3),
];

impl<T, I, C> GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Debug + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
{
    /// Sample and apply a one-qubit Pauli channel.
    ///
    /// Return codes are `0 = I`, `1 = X`, `2 = Y`, `3 = Z`.
    pub fn sample_apply_pauli_error_choice(&mut self, addr: usize, p: [f64; 3]) -> u8 {
        debug_assert!(p.iter().all(|p| p.is_finite() && (0.0..=1.0).contains(p)));
        let choice = self.sample_choice(&p);
        self.apply_pauli_choice(addr, choice);
        choice
    }

    /// Sample and apply a one-qubit depolarizing channel.
    ///
    /// Return codes are `0 = I`, `1 = X`, `2 = Y`, `3 = Z`.
    pub fn sample_apply_depolarize1_choice(&mut self, addr: usize, p: f64) -> u8 {
        self.sample_apply_pauli_error_choice(addr, [p / 3.0, p / 3.0, p / 3.0])
    }

    /// Sample and apply a two-qubit Pauli channel.
    ///
    /// Return codes are the Pauli-pair index where `0 = II`, `1 = IX`,
    /// `2 = IY`, ..., `15 = ZZ`. If either qubit is already lost, this matches
    /// the normal channel implementation by consuming no RNG and returning
    /// `0`.
    pub fn sample_apply_two_qubit_pauli_error_choice(
        &mut self,
        addr0: usize,
        addr1: usize,
        p: [f64; 15],
    ) -> u8 {
        if self.is_lost[addr0] || self.is_lost[addr1] {
            return 0;
        }
        debug_assert!(p.iter().all(|p| p.is_finite() && (0.0..=1.0).contains(p)));
        let choice = self.sample_choice(&p);
        if choice != 0 {
            self.apply_two_qubit_pauli_choice(addr0, addr1, choice);
        }
        choice
    }

    /// Sample and apply a two-qubit depolarizing channel.
    ///
    /// Return codes are the same as
    /// [`sample_apply_two_qubit_pauli_error_choice`](Self::sample_apply_two_qubit_pauli_error_choice).
    pub fn sample_apply_depolarize2_choice(&mut self, addr0: usize, addr1: usize, p: f64) -> u8 {
        self.sample_apply_two_qubit_pauli_error_choice(addr0, addr1, [p / 15.0; 15])
    }

    /// Sample and apply a single-qubit loss channel.
    ///
    /// Return codes are variable-length:
    /// - `[0]` means the qubit survived.
    /// - `[1, m]` means it was lost after collapse/reset, where `m` is
    ///   `0 = collapsed to |0>`, `1 = collapsed to |1>`, `2 = already lost`.
    pub fn sample_apply_loss_choice(&mut self, addr: usize, p: f64) -> Vec<u8> {
        debug_assert!(p.is_finite() && (0.0..=1.0).contains(&p));
        if p < self.tableau.rng.random::<f64>() {
            return vec![0];
        }

        let collapse = self.sample_apply_reset_choice(addr);
        self.is_lost[addr] = true;
        vec![1, collapse]
    }

    /// Sample and apply a correlated two-qubit loss channel.
    ///
    /// Return codes are variable-length:
    /// - `[0]`: no loss.
    /// - `[1, m0, m1]`: both qubits lost; `m*` are reset-collapse codes.
    /// - `[2, m0]`: only `addr0` lost.
    /// - `[3, m1]`: only `addr1` lost.
    /// - `[4, ...]`: `addr0` was already lost, followed by a single-loss code
    ///   for `addr1` with probability `p[2]`.
    /// - `[5, ...]`: `addr1` was already lost, followed by a single-loss code
    ///   for `addr0` with probability `p[2]`.
    pub fn sample_apply_correlated_loss_choice(
        &mut self,
        addr0: usize,
        addr1: usize,
        p: [f64; 3],
    ) -> Vec<u8> {
        debug_assert!(p.iter().all(|p| p.is_finite() && (0.0..=1.0).contains(p)));
        if self.is_lost[addr0] {
            let mut choice = vec![4];
            choice.extend(self.sample_apply_loss_choice(addr1, p[2]));
            return choice;
        }
        if self.is_lost[addr1] {
            let mut choice = vec![5];
            choice.extend(self.sample_apply_loss_choice(addr0, p[2]));
            return choice;
        }

        let r = self.tableau.rng.random::<f64>();
        if r < p[0] {
            let m0 = self.sample_apply_reset_choice(addr0);
            let m1 = self.sample_apply_reset_choice(addr1);
            self.is_lost[addr0] = true;
            self.is_lost[addr1] = true;
            return vec![1, m0, m1];
        }
        if r < p[0] + p[1] {
            if self.tableau.rng.random::<bool>() {
                let m1 = self.sample_apply_reset_choice(addr1);
                self.is_lost[addr1] = true;
                vec![3, m1]
            } else {
                let m0 = self.sample_apply_reset_choice(addr0);
                self.is_lost[addr0] = true;
                vec![2, m0]
            }
        } else {
            vec![0]
        }
    }

    fn sample_choice(&mut self, probabilities: &[f64]) -> u8 {
        let r = self.tableau.rng.random::<f64>();
        let mut cumulative = 0.0;
        for (i, p) in probabilities.iter().enumerate() {
            cumulative += *p;
            if cumulative > r {
                return (i + 1) as u8;
            }
        }
        0
    }

    fn apply_pauli_choice(&mut self, addr: usize, choice: u8) {
        match choice {
            0 => {}
            1 => self.x(addr),
            2 => self.y(addr),
            3 => self.z(addr),
            _ => unreachable!("single-qubit Pauli choice must be in 0..=3"),
        }
    }

    fn apply_two_qubit_pauli_choice(&mut self, addr0: usize, addr1: usize, choice: u8) {
        let (p0, p1) = PAULI_PAIRS[choice as usize];
        self.apply_pauli_choice(addr0, p0);
        self.apply_pauli_choice(addr1, p1);
    }

    fn sample_apply_reset_choice(&mut self, addr: usize) -> u8 {
        let outcome = self.measure(addr);
        self.measurement_record.pop();
        if let Some(true) = outcome {
            self.x(addr);
        }
        measurement_code(outcome)
    }
}

fn measurement_code(outcome: Option<bool>) -> u8 {
    match outcome {
        Some(false) => 0,
        Some(true) => 1,
        None => 2,
    }
}
