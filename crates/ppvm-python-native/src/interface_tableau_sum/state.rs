// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[rustfmt::skip]
macro_rules! create_sum_state {
    ($tab_name: ident, $sampler_name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pyclass]
        pub struct $tab_name {
            inner: $type,
        }
        #[pymethods]
        impl $tab_name {
            #[new]
            #[pyo3(signature = (
                                                                        n_qubits,
                                                                        min_abs_coeff = 1e-10,
                                                                        sum_cutoff = 1e-8,
                                                                        seed = None,
                                                                    ))]
            pub fn new(
                n_qubits: usize,
                min_abs_coeff: f64,
                sum_cutoff: f64,
                seed: Option<u64>,
            ) -> Self {
                let inner: $type = match seed {
                    Some(s) => {
                        <$type>::new_with_seed(n_qubits, min_abs_coeff, sum_cutoff, s)
                    }
                    None => <$type>::new(n_qubits, min_abs_coeff, sum_cutoff),
                };
                Self { inner }
            }

            /// Number of branches currently in the sum.
            pub fn __len__(&self) -> usize {
                self.inner.len()
            }

            pub fn len(&self) -> usize {
                self.inner.len()
            }

            pub fn is_empty(&self) -> bool {
                self.inner.is_empty()
            }

            /// Mid-circuit measurement probabilities `(p_zero, p_one, p_lost)`.
            ///
            /// Unlike `GeneralizedTableau.measure` (which samples and returns a
            /// single outcome), the sum-form measurement branches each entry
            /// into its three Z-basis outcomes and returns the aggregated
            /// probabilities — use this for analytic measurement statistics,
            /// and a `Sampler` from `.sampler()` for stochastic samples.
            pub fn measure(&mut self, addr0: usize) -> (f64, f64, f64) {
                self.inner.measure(addr0)
            }
            /// Force a truncation pass.
            ///
            /// Truncation is already applied automatically inside every noise
            /// op; this is exposed for users who want to prune after a long
            /// run of unitary gates.
            pub fn truncate(&mut self) {
                self.inner.truncate();
            }

            /// Compile a `Sampler` snapshotting the current state.
            ///
            /// The returned sampler holds its own RNG and a sorted copy of the
            /// sum's branches; further gates / noise on this tableau do not
            /// affect it. Call `.sample_shots(n)` on the result to draw shots.
            /// Two samplers compiled in a row use independent RNG sequences.
            pub fn sampler(&mut self) -> $sampler_name {
                $sampler_name {
                    inner: self.inner.sampler(),
                }
            }

            fn __copy__(&self) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }

            fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }
        }
    };
}
