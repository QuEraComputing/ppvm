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
    pub(super) fn draw(&mut self, probabilities: &[f64]) -> Option<usize> {
        debug_assert!(probabilities.iter().all(|&p| p >= 0.0));
        debug_assert!(probabilities.iter().sum::<f64>() <= 1.0 + 1e-9);
        let r = self.rng.random::<f64>();
        let mut cumulative = 0.0;
        probabilities.iter().position(|&p| {
            cumulative += p;
            r < cumulative
        })
    }

    /// Draw event `k` with probability `probabilities[k] / total`, given that
    /// one happens. Rounding falls back to the last possible event.
    fn draw_given_event(&mut self, probabilities: &[f64], total: f64) -> usize {
        let r = self.rng.random::<f64>() * total;
        let mut cumulative = 0.0;
        probabilities
            .iter()
            .position(|&p| {
                cumulative += p;
                r < cumulative
            })
            .unwrap_or_else(|| probabilities.iter().rposition(|&p| p > 0.0).unwrap())
    }

    /// First index `≥ i` that gets an error when each target does independently
    /// with probability `p`: a geometric skip from one random number.
    fn next_error(&mut self, i: usize, p: f64) -> usize {
        if p <= 0.0 {
            return usize::MAX;
        }
        let u: f64 = self.rng.random();
        // ln(1 − u) / ln(1 − p) ≥ 0; the cast saturates, so huge skips end the loop.
        let skip = ((-u).ln_1p() / (-p).ln_1p()).floor() as usize;
        i.saturating_add(skip)
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

/// `X`, `Y`, `Z` with probabilities `p[0]`, `p[1]`, `p[2]`. The batched
/// methods skip straight to the next error instead of drawing per target.
impl<V: VirtualState> PauliError<f64> for GeneralizedTableau<V> {
    fn pauli_error(&mut self, q: usize, p: [f64; 3]) {
        if let Some(k) = self.draw(&p) {
            self.pauli(q, k + 1);
        }
    }

    fn pauli_error_many(&mut self, targets: &[usize], p: [f64; 3]) {
        let total: f64 = p.iter().sum();
        let mut i = self.next_error(0, total);
        while i < targets.len() {
            let k = self.draw_given_event(&p, total);
            self.pauli(targets[i], k + 1);
            i = self.next_error(i + 1, total);
        }
    }

    fn x_error_many(&mut self, targets: &[usize], p: f64) {
        self.pauli_error_many(targets, [p, 0.0, 0.0]);
    }

    fn y_error_many(&mut self, targets: &[usize], p: f64) {
        self.pauli_error_many(targets, [0.0, p, 0.0]);
    }

    fn z_error_many(&mut self, targets: &[usize], p: f64) {
        self.pauli_error_many(targets, [0.0, 0.0, p]);
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

    /// Skip sampling over pairs, as for `pauli_error_many`.
    fn two_qubit_pauli_error_many(&mut self, pairs: &[(usize, usize)], p: [f64; 15]) {
        let total: f64 = p.iter().sum();
        let mut i = self.next_error(0, total);
        while i < pairs.len() {
            let k = self.draw_given_event(&p, total);
            let (a, b) = pairs[i];
            self.pauli(a, (k + 1) / 4);
            self.pauli(b, (k + 1) % 4);
            i = self.next_error(i + 1, total);
        }
    }
}

/// Stim's `DEPOLARIZE1(p)`: each of `X`, `Y`, `Z` with probability `p/3`.
impl<V: VirtualState> Depolarizing<f64> for GeneralizedTableau<V> {
    fn depolarize1(&mut self, q: usize, p: f64) {
        self.pauli_error(q, [p / 3.0; 3]);
    }

    fn depolarize1_many(&mut self, targets: &[usize], p: f64) {
        self.pauli_error_many(targets, [p / 3.0; 3]);
    }
}

/// Stim's `DEPOLARIZE2(p)`: each of the 15 non-identity Paulis with
/// probability `p/15`.
impl<V: VirtualState> Depolarizing2<f64> for GeneralizedTableau<V> {
    fn depolarize2(&mut self, a: usize, b: usize, p: f64) {
        self.two_qubit_pauli_error(a, b, [p / 15.0; 15]);
    }

    fn depolarize2_many(&mut self, pairs: &[(usize, usize)], p: f64) {
        self.two_qubit_pauli_error_many(pairs, [p / 15.0; 15]);
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

    /// Batched channels act on each target independently with the channel's
    /// rates: per Pauli, at either end of the batch, and on neighbouring pairs.
    #[test]
    fn batched_channels_match_independent_rates() {
        let (n, shots, p) = (400, 100, 0.06);
        let targets: Vec<usize> = (0..n).collect();
        let samples = |count: usize| (count * shots) as f64;
        let within = |hits: usize, trials: f64, rate: f64| {
            let sigma = (rate * (1.0 - rate) / trials).sqrt();
            ((hits as f64 / trials) - rate).abs() <= 5.0 * sigma
        };
        let channels: [(&str, [f64; 3], fn(&mut Sim, &[usize], f64)); 4] = [
            ("X", [p, 0.0, 0.0], |s, t, p| s.x_error_many(t, p)),
            ("Y", [0.0, p, 0.0], |s, t, p| s.y_error_many(t, p)),
            ("Z", [0.0, 0.0, p], |s, t, p| s.z_error_many(t, p)),
            ("DEPOLARIZE1", [p / 3.0; 3], |s, t, p| {
                s.depolarize1_many(t, p)
            }),
        ];
        for (name, rates, channel) in channels {
            let (mut per_pauli, mut ends, mut pairs) = ([0; 3], [0; 2], 0);
            for seed in 0..shots {
                let got = applied(n, seed as u64, |s| channel(s, &targets, p));
                let hit: Vec<bool> = got.iter().map(|&b| b != bits('I')).collect();
                for b in &got {
                    if let Some(k) = "XYZ".chars().position(|c| bits(c) == *b) {
                        per_pauli[k] += 1;
                    }
                }
                ends[0] += hit[..n / 10].iter().filter(|&&h| h).count();
                ends[1] += hit[n - n / 10..].iter().filter(|&&h| h).count();
                pairs += hit.windows(2).filter(|w| w[0] && w[1]).count();
            }
            for k in 0..3 {
                assert!(
                    within(per_pauli[k], samples(n), rates[k]),
                    "{name}: {per_pauli:?}"
                );
            }
            for end in ends {
                assert!(within(end, samples(n / 10), p), "{name}: ends {ends:?}");
            }
            assert!(
                within(pairs, samples(n - 1), p * p),
                "{name}: pairs {pairs}"
            );
        }
    }

    /// `p = 0` never applies the Pauli, `p = 1` applies it to every target.
    #[test]
    fn batched_channels_at_zero_and_one() {
        let targets: Vec<usize> = (0..50).collect();
        let none = applied(50, 0, |s| s.x_error_many(&targets, 0.0));
        assert!(none.iter().all(|&b| b == bits('I')));
        let all = applied(50, 0, |s| s.x_error_many(&targets, 1.0));
        assert!(all.iter().all(|&b| b == bits('X')));
    }

    /// Batched two-qubit channels act on each pair independently: each of the
    /// 15 Paulis at its rate (all distinct, so a mixed-up order shows), at
    /// either end of the batch, and on neighbouring pairs.
    #[test]
    fn batched_two_qubit_channels_match_independent_rates() {
        let (n, shots) = (400, 100);
        let pairs: Vec<(usize, usize)> = (0..n / 2).map(|i| (2 * i, 2 * i + 1)).collect();
        let rates: [f64; 15] = std::array::from_fn(|k| 0.0005 * (k + 1) as f64);
        let total: f64 = rates.iter().sum();
        let order: Vec<(bool, bool, bool, bool)> = "IX IY IZ XI XX XY XZ YI YX YY YZ ZI ZX ZY ZZ"
            .split(' ')
            .map(|pair| {
                let c: Vec<char> = pair.chars().collect();
                let ((xa, za), (xb, zb)) = (bits(c[0]), bits(c[1]));
                (xa, za, xb, zb)
            })
            .collect();
        let samples = |count: usize| (count * shots) as f64;
        let within = |hits: usize, trials: f64, rate: f64| {
            let sigma = (rate * (1.0 - rate) / trials).sqrt();
            ((hits as f64 / trials) - rate).abs() <= 5.0 * sigma
        };
        let channels: [(&str, [f64; 15], f64, fn(&mut Sim, &[(usize, usize)])); 2] = [
            ("PAULI_CHANNEL_2", rates, total, |s, t| {
                s.two_qubit_pauli_error_many(t, std::array::from_fn(|k| 0.0005 * (k + 1) as f64))
            }),
            ("DEPOLARIZE2", [0.06 / 15.0; 15], 0.06, |s, t| {
                s.depolarize2_many(t, 0.06)
            }),
        ];
        for (name, rates, total, channel) in channels {
            let (mut per_pauli, mut ends, mut both) = ([0; 15], [0; 2], 0);
            let tenth = pairs.len() / 10;
            for seed in 0..shots {
                let got = applied(n, seed as u64, |s| channel(s, &pairs));
                let kinds: Vec<Option<usize>> = got
                    .chunks(2)
                    .map(|q| {
                        order
                            .iter()
                            .position(|&o| o == (q[0].0, q[0].1, q[1].0, q[1].1))
                    })
                    .collect();
                for k in kinds.iter().flatten() {
                    per_pauli[*k] += 1;
                }
                let hit: Vec<bool> = kinds.iter().map(Option::is_some).collect();
                ends[0] += hit[..tenth].iter().filter(|&&h| h).count();
                ends[1] += hit[hit.len() - tenth..].iter().filter(|&&h| h).count();
                both += hit.windows(2).filter(|w| w[0] && w[1]).count();
            }
            for k in 0..15 {
                assert!(
                    within(per_pauli[k], samples(pairs.len()), rates[k]),
                    "{name}: {per_pauli:?}"
                );
            }
            for end in ends {
                assert!(within(end, samples(tenth), total), "{name}: ends {ends:?}");
            }
            assert!(
                within(both, samples(pairs.len() - 1), total * total),
                "{name}: both {both}"
            );
        }
    }
}
