// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use rand::RngExt;

use ppvm_traits_2::Measure;

use super::GeneralizedTableau;
use crate::VirtualState;

impl<V: VirtualState> GeneralizedTableau<V> {
    /// Unnormalized weights of the outcomes `0` and `1` of `Z_q`.
    pub(super) fn z_weights(&mut self, q: usize) -> [f64; 2] {
        let mut view = self.tableau.virtual_view();
        self.state.z_weights(&mut view, q)
    }

    /// Project onto outcome `m` of `Z_q`, which must have nonzero weight.
    pub(super) fn collapse(&mut self, q: usize, m: bool) {
        let mut view = self.tableau.virtual_view();
        self.state.collapse(&mut view, q, m);
    }

    /// Sample `Z_q` and collapse, whatever the status of `q`.
    pub(super) fn sample_z(&mut self, q: usize) -> bool {
        // Hot path: a stabilizer state measures on the tableau alone.
        if self.state.is_zero() {
            let rng = &mut self.rng;
            return self.tableau.measure(q, || rng.random());
        }
        let [w0, w1] = self.z_weights(q);
        let m = self.rng.random::<f64>() * (w0 + w1) < w1;
        self.collapse(q, m);
        m
    }
}

/// Measure `Z_q` and collapse; `Some(true)` for `|1⟩`. A lost qubit gives
/// `None`, a leaked one its pinned bit.
impl<V: VirtualState> Measure for GeneralizedTableau<V> {
    fn measure(&mut self, q: usize) -> Option<bool> {
        if self.is_lost(q) {
            return None;
        }
        Some(self.sample_z(q))
    }
}

#[cfg(test)]
mod tests {
    use ppvm_traits_2::{Clifford, Reset};

    use super::*;
    use crate::{ActiveDormant, DenseState};

    type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

    /// Bell pair: the two outcomes always agree, and both occur.
    #[test]
    fn bell_pair_outcomes_agree() {
        let mut ones = 0;
        for seed in 0..200 {
            let mut s = Sim::new_with_seed(2, seed);
            s.h(0);
            s.cnot(0, 1);
            let (a, b) = (s.measure(0).unwrap(), s.measure(1).unwrap());
            assert_eq!(a, b);
            ones += usize::from(a);
        }
        assert!((60..140).contains(&ones), "{ones} of 200");
    }

    /// A measured qubit stays in its outcome, and reset leaves `|0⟩`.
    #[test]
    fn measurement_repeats_and_reset_clears() {
        for seed in 0..50 {
            let mut s = Sim::new_with_seed(3, seed);
            s.h(0);
            s.cnot(0, 1);
            s.sqrt_x(2);
            s.cz(1, 2);
            for q in 0..3 {
                let m = s.measure(q);
                assert_eq!(s.measure(q), m);
            }
            for q in 0..3 {
                s.reset(q);
                assert_eq!(s.measure(q), Some(false));
            }
        }
    }

    /// A fork carries the collapsed state: after measuring one half of a Bell
    /// pair, the fork's other half agrees with it.
    #[test]
    fn fork_keeps_the_collapsed_state() {
        for seed in 0..50 {
            let mut s = Sim::new_with_seed(2, seed);
            s.h(0);
            s.cnot(0, 1);
            let a = s.measure(0);
            let mut f = s.fork(Some(seed + 1000));
            assert_eq!(f.measure(1), a);
        }
    }

    /// Forks of a superposition sample independently: one seed repeats its
    /// outcome, different seeds give both outcomes.
    #[test]
    fn forks_sample_independently() {
        let mut s = Sim::new_with_seed(1, 0);
        s.h(0);
        let mut ones = 0;
        for seed in 0..200 {
            let m = s.fork(Some(seed)).measure(0);
            assert_eq!(s.fork(Some(seed)).measure(0), m);
            ones += usize::from(m.unwrap());
        }
        assert!((60..140).contains(&ones), "{ones} of 200");
    }

    /// `X` then measure is deterministic: `Z|1⟩ = −|1⟩`.
    #[test]
    fn x_flips_the_outcome() {
        let mut s = Sim::new_with_seed(1, 0);
        s.x(0);
        assert_eq!(s.measure(0), Some(true));
    }
}
