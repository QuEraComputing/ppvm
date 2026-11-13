// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Pauli channels as trajectories: draw one Pauli and apply it as a physical
//! gate, so `|χ⟩` is untouched.

use rand::RngExt;

use ppvm_traits_2::{Clifford, Depolarizing, Depolarizing2, PauliError, TwoQubitPauliError};

use super::GeneralizedTableau;
use crate::VirtualState;

impl<V: VirtualState> GeneralizedTableau<V> {
    /// Draw event `k` with probability `probabilities[k]`, or `None` with the
    /// remaining probability. Uses one random number.
    fn draw(&mut self, probabilities: &[f64]) -> Option<usize> {
        debug_assert!(probabilities.iter().all(|&p| p >= 0.0));
        debug_assert!(probabilities.iter().sum::<f64>() <= 1.0 + 1e-9);
        let r = self.rng.random::<f64>();
        let mut cumulative = 0.0;
        probabilities.iter().position(|&p| {
            cumulative += p;
            r < cumulative
        })
    }

    /// Apply Pauli `k` (0 = I, 1 = X, 2 = Y, 3 = Z) to `q`.
    fn pauli(&mut self, q: usize, k: usize) {
        match k {
            0 => {}
            1 => self.x(q),
            2 => self.y(q),
            _ => self.z(q),
        }
    }
}

/// `X`, `Y`, `Z` with probabilities `p[0]`, `p[1]`, `p[2]`.
impl<V: VirtualState> PauliError<f64> for GeneralizedTableau<V> {
    fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        if let Some(k) = self.draw(&p) {
            self.pauli(q, k + 1);
        }
    }
}

/// `p[k]` is the probability of `P_a ⊗ P_b` with `4a + b = k + 1`, in the
/// trait's order `IX, IY, IZ, XI, …, ZZ`.
impl<V: VirtualState> TwoQubitPauliError<f64> for GeneralizedTableau<V> {
    fn two_qubit_pauli_error(&mut self, a: usize, b: usize, p: [f64; 15]) {
        if let Some(k) = self.draw(&p) {
            self.pauli(a, (k + 1) / 4);
            self.pauli(b, (k + 1) % 4);
        }
    }
}

/// Stim's `DEPOLARIZE1(p)`: each of `X`, `Y`, `Z` with probability `p/3`.
impl<V: VirtualState> Depolarizing<f64> for GeneralizedTableau<V> {
    fn depolarize1(&mut self, q: usize, p: f64) {
        self.pauli_error(q, [p / 3.0; 3]);
    }
}

/// Stim's `DEPOLARIZE2(p)`: each of the 15 non-identity Paulis with
/// probability `p/15`.
impl<V: VirtualState> Depolarizing2<f64> for GeneralizedTableau<V> {
    fn depolarize2(&mut self, a: usize, b: usize, p: f64) {
        self.two_qubit_pauli_error(a, b, [p / 15.0; 15]);
    }
}

#[cfg(test)]
mod tests {
    use ppvm_traits_2::Measure;

    use super::*;
    use crate::{ActiveDormant, DenseState};

    type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

    /// The `(x, z)` bits of the Pauli applied to each qubit of a fresh `n`-qubit
    /// state by `channel`: X bits flip `|0⟩`, Z bits flip `|+⟩`.
    fn applied(n: usize, seed: u64, channel: impl Fn(&mut Sim)) -> Vec<(bool, bool)> {
        let (mut zero, mut plus) = (Sim::new_with_seed(n, seed), Sim::new_with_seed(n, seed));
        (0..n).for_each(|q| plus.h(q));
        channel(&mut zero);
        channel(&mut plus);
        (0..n)
            .map(|q| {
                plus.h(q);
                (zero.measure(q).unwrap(), plus.measure(q).unwrap())
            })
            .collect()
    }

    fn bits(name: char) -> (bool, bool) {
        match name {
            'I' => (false, false),
            'X' => (true, false),
            'Y' => (true, true),
            _ => (false, true),
        }
    }

    /// With all probability on one entry, exactly that Pauli is applied, in
    /// the order the traits document.
    #[test]
    fn one_hot_channels_apply_the_documented_pauli() {
        for (k, name) in "XYZ".chars().enumerate() {
            let mut p = [0.0; 3];
            p[k] = 1.0;
            assert_eq!(applied(1, 0, |s| s.pauli_error(0, p)), [bits(name)]);
        }
        let order = "IX IY IZ XI XX XY XZ YI YX YY YZ ZI ZX ZY ZZ";
        for (k, pair) in order.split(' ').enumerate() {
            let mut p = [0.0; 15];
            p[k] = 1.0;
            let names: Vec<char> = pair.chars().collect();
            assert_eq!(
                applied(2, 0, |s| s.two_qubit_pauli_error(0, 1, p)),
                [bits(names[0]), bits(names[1])],
                "{pair}"
            );
        }
    }

    /// Flip frequencies on `|0…0⟩` match the channels: `DEPOLARIZE1(p)` flips
    /// with `2p/3` (X, Y), `DEPOLARIZE2(p)` flips a given qubit with `8p/15`
    /// and both with `4p/15`.
    #[test]
    fn depolarizing_flip_rates() {
        let (p, shots) = (0.3, 20_000);
        let within = |count: usize, rate: f64| {
            let sigma = (rate * (1.0 - rate) / shots as f64).sqrt();
            ((count as f64 / shots as f64) - rate).abs() < 5.0 * sigma
        };
        let mut flips = 0;
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(1, seed as u64);
            s.depolarize1(0, p);
            flips += usize::from(s.measure(0).unwrap());
        }
        assert!(within(flips, 2.0 * p / 3.0), "{flips} of {shots}");
        let (mut first, mut both) = (0, 0);
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(2, seed as u64);
            s.depolarize2(0, 1, p);
            let (a, b) = (s.measure(0).unwrap(), s.measure(1).unwrap());
            first += usize::from(a);
            both += usize::from(a && b);
        }
        assert!(within(first, 8.0 * p / 15.0), "{first} of {shots}");
        assert!(within(both, 4.0 * p / 15.0), "{both} of {shots}");
    }
}
