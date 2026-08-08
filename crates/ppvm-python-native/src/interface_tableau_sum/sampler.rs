// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_sum_sampler {
    ($tab_name: ident, $sampler_name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pyclass]
        pub struct $sampler_name {
            inner: crate::backend::MixtureSampler<{ $storage }, $indexType>,
        }

        #[pymethods]
        impl $sampler_name {
            /// Draw a single full-register sample. Per-qubit outcomes are
            /// encoded as `0 = |0>`, `1 = |1>`, `2 = lost`.
            pub fn sample(&mut self) -> Vec<u8> {
                self.inner
                    .sample()
                    .into_iter()
                    .map(measurement_to_u8)
                    .collect()
            }

            /// Draw `num_shots` full-register samples in parallel.
            ///
            /// Runs on the rayon thread pool with the GIL released, so it
            /// scales across cores for batched sampling.
            pub fn sample_shots(&mut self, py: Python<'_>, num_shots: usize) -> Vec<Vec<u8>> {
                let raw = py.detach(|| self.inner.sample_shots(num_shots));
                raw.into_iter()
                    .map(|shot| shot.into_iter().map(measurement_to_u8).collect())
                    .collect()
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
