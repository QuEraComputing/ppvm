// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Allocating predictor-corrector step, built from the map-based leakage and
//! the column-assembled exponential, as the numerical reference for the
//! buffer-reusing [`LindbladSpec::pc_step`].

use crate::truncate::{add_leakage_capped, cap_basis, prune_basis};
use crate::{Error, LindbladSpec, PcStepConfig, Word, mf_expm};

impl<const C: usize> LindbladSpec<C> {
    pub(crate) fn pc_step_reference(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        self.run_in_pool(cfg, |this| {
            this.pc_step_reference_inner(basis, coeffs, dt, protected, cfg)
        })
    }

    fn pc_step_reference_inner(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        let PcStepConfig {
            max_basis,
            admit_basis,
            drop_tol,
            tau_add,
            ..
        } = *cfg;
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);

        let leak = self.leakage_with_prune(basis, coeffs, protected, admit, tau_add)?;
        add_leakage_capped(basis, coeffs, leak, admit);
        let coeffs_predict = mf_expm::expm_apply_mf_reference(self, basis, dt, coeffs, drop_tol);

        let n_predict = basis.len();
        if admit > n_predict {
            let leak2 =
                self.leakage_with_prune(basis, &coeffs_predict, protected, admit, tau_add)?;
            add_leakage_capped(basis, coeffs, leak2, admit);
        }
        if basis.len() == n_predict {
            *coeffs = coeffs_predict;
        } else {
            drop(coeffs_predict);
            *coeffs = mf_expm::expm_apply_mf_reference(self, basis, dt, coeffs, drop_tol);
        }

        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(())
    }
}
