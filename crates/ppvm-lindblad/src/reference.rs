// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Allocation-independent predictor/corrector reference for numerical regression tests.

use crate::step::Phase;
use crate::truncate::{
    add_leakage_capped, cap_basis, cap_map_to_room, order_by_desc_mag, prune_basis,
};
use crate::word::word_hash;
use crate::{Error, LindbladSpec, PcStepConfig, PcStepTimings, Word};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use rayon::prelude::*;
const CHUNK_SIZE: usize = 4096;

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
            this.pc_step_reference_inner(basis, coeffs, dt, protected, cfg, false)
                .map(|_| ())
        })
    }
    fn pc_step_reference_inner(
        &self,
        basis: &mut Vec<Word<C>>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word<C>],
        cfg: &PcStepConfig,
        timed: bool,
    ) -> Result<PcStepTimings, Error> {
        let PcStepConfig {
            max_basis,
            admit_basis,
            drop_tol,
            tau_add,
            ..
        } = *cfg;
        // Admission bound: enrichment may grow the live basis to `admit`
        // >= `max_basis`; the final `cap_basis` then keeps the top-
        // `max_basis` strings by evolved |coeff| over the whole union
        // (retained + admitted) — rank displacement. With `admit_basis =
        // None` admission is bounded by `max_basis` itself, `cap_basis` is
        // a no-op, and membership turnover requires `drop_tol > 0`.
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);
        let mut t = PcStepTimings::default();

        // 1. First-hop expansion. After this, `coeffs` contains the pre-step
        // coefficients followed by zeros for the newly-added leakage strings.
        // We rely on `coeffs` itself as the pre-step buffer for the corrector
        // — the reference exponential only borrows its input coefficients.
        let p = Phase::start(timed);
        let leak = self.leakage_reference_with_prune(basis, coeffs, protected, admit, tau_add)?;
        p.stop(&mut t.leakage1_us);

        let p = Phase::start(timed);
        add_leakage_capped(basis, coeffs, leak, admit);
        p.stop(&mut t.expand1_us);

        // 2. Predictor: the reference exponential reads `coeffs` and returns a
        // new owned vector with the predicted state.
        let p = Phase::start(timed);
        let coeffs_predict =
            crate::mf_expm::expm_apply_mf_reference(self, basis, dt, coeffs, drop_tol);
        p.stop(&mut t.expm1_us);

        // 3. Second-hop expansion from the predicted state. Extend `coeffs`
        // with zeros for any newly-added second-hop strings so it remains a
        // valid input (pre-step state) for the corrector. Once the basis is
        // full, the first hop usually fills the whole admission room; with
        // `room = 0` the second hop can admit nothing, so its leakage pass
        // is skipped.
        let n_predict = basis.len();
        if admit > n_predict {
            let p = Phase::start(timed);
            let leak2 = self.leakage_reference_with_prune(
                basis,
                &coeffs_predict,
                protected,
                admit,
                tau_add,
            )?;
            p.stop(&mut t.leakage2_us);

            let p = Phase::start(timed);
            add_leakage_capped(basis, coeffs, leak2, admit);
            p.stop(&mut t.expand2_us);
        }

        // 4. Corrector: redo from pre-step state on the doubly-enlarged basis.
        // If the second hop admitted nothing, the corrector would repeat the
        // predictor's computation exactly, so the predicted state is kept.
        if basis.len() == n_predict {
            *coeffs = coeffs_predict;
        } else {
            drop(coeffs_predict);
            let p = Phase::start(timed);
            *coeffs = crate::mf_expm::expm_apply_mf_reference(self, basis, dt, coeffs, drop_tol);
            p.stop(&mut t.expm2_us);
        }

        // 5. Prune basis entries below `drop_tol` (protected words never dropped).
        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(t)
    }

    pub(crate) fn leakage_reference_with_prune(
        &self,
        basis: &[Word<C>],
        coeffs: &[f64],
        protected: &[Word<C>],
        max_basis: usize,
        tau_add: f64,
    ) -> Result<Vec<(Word<C>, f64)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        // Hash-only membership tables: storing 8-byte `u64` keys instead
        // of 48-byte Words shrinks the in-basis structure ~6×, keeping it
        // in L3 (and often L2) at basis sizes where the full-Word version
        // would spill to DRAM.
        let in_basis: FxHashMap<u64, ()> = basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        let order = order_by_desc_mag(coeffs);
        let room = max_basis.saturating_sub(basis.len());
        let n_qubits = self.n_qubits();
        let mut merged: FxHashMap<Word<C>, f64> = FxHashMap::default();
        for chunk_indices in order.chunks(CHUNK_SIZE) {
            let local: Vec<Vec<(Word<C>, f64)>> = chunk_indices
                .par_iter()
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word<C>, Complex<f64>>::with_capacity_and_hasher(
                                128,
                                FxBuildHasher::default(),
                            ),
                        )
                    },
                    |(s1, s2, lm), &i| {
                        let p = &basis[i];
                        let c = coeffs[i];
                        let terms = self.compute_action_terms(p, s1, s2, lm);
                        let mut out = Vec::with_capacity(terms.len());
                        for (w, v) in terms.iter() {
                            let h = word_hash(w);
                            if !in_basis.contains_key(&h) && !protected_set.contains_key(&h) {
                                out.push((*w, c * *v));
                            }
                        }
                        out
                    },
                )
                .collect();
            for v in local {
                for (k, val) in v {
                    *merged.entry(k).or_insert(0.0) += val;
                }
            }
            cap_map_to_room(&mut merged, room);
        }
        // Rate-based admission: keep only candidates whose leakage rate
        // exceeds `tau_add`. `tau_add = 0` admits everything except exact
        // zeros.
        Ok(merged
            .into_iter()
            .filter(|(_, c)| c.abs() > tau_add)
            .collect())
    }
}
