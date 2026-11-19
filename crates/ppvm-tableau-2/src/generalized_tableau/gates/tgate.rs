// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num_complex::Complex64;
use ppvm_traits_2::TGate;

use super::super::GeneralizedTableau;
use crate::VirtualState;

/// `T = diag(1, ω) = ½(1 + ω)·I + ½(1 − ω)·Z` with `ω = e^{iπ/4}`, exact
/// including the global phase.
impl<V: VirtualState> TGate for GeneralizedTableau<V> {
    fn t(&mut self, q: usize) {
        let w = Complex64::from_polar(1.0, std::f64::consts::FRAC_PI_4);
        self.rotate_z(&[q], (1.0 + w) / 2.0, (1.0 - w) / 2.0);
    }

    fn t_dag(&mut self, q: usize) {
        let w = Complex64::from_polar(1.0, -std::f64::consts::FRAC_PI_4);
        self.rotate_z(&[q], (1.0 + w) / 2.0, (1.0 - w) / 2.0);
    }
}
