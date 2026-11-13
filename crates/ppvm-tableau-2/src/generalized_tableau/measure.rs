// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use rand::RngExt;

use ppvm_traits_2::Measure;

use super::GeneralizedTableau;
use crate::VirtualState;

/// Measure `Z_q` and collapse; `Some(true)` for `|1⟩`. Every qubit is live
/// until loss is modelled, so the outcome is never `None`.
impl<V: VirtualState> Measure for GeneralizedTableau<V> {
    fn measure(&mut self, q: usize) -> Option<bool> {
        let mut view = self.tableau.virtual_view();
        let [w0, w1] = self.state.z_weights(&mut view, q);
        let m = self.rng.random::<f64>() * (w0 + w1) < w1;
        self.state.collapse(&mut view, q, m);
        Some(m)
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

    /// `X` then measure is deterministic: `Z|1⟩ = −|1⟩`.
    #[test]
    fn x_flips_the_outcome() {
        let mut s = Sim::new_with_seed(1, 0);
        s.x(0);
        assert_eq!(s.measure(0), Some(true));
    }
}
