// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_tableau_stim {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pymethods]
        impl $name {
            pub fn run(
                &mut self,
                prog: &crate::stim_program::PyStimProgram,
            ) -> pyo3::PyResult<Vec<u8>> {
                let mut results = Vec::with_capacity(prog.measurement_count());
                ppvm_stim::execute_validated(&prog.instructions, &mut self.inner, &mut results);
                Ok(results
                    .into_iter()
                    .map(crate::interface_tableau::measurement_to_u8)
                    .collect())
            }

            /// Multi-shot sampling: builds a fresh tableau per shot.
            ///
            /// Shots run in parallel on rayon's global thread pool (GIL
            /// released), falling back to serial for small batches. Shot `i`
            /// is seeded with `seed.wrapping_add(i)` when `seed` is given
            /// (wrapping mod 2⁶⁴), so results are reproducible and
            /// independent of the thread count; set the `RAYON_NUM_THREADS`
            /// environment variable to control the pool size.
            #[staticmethod]
            #[pyo3(signature = (prog, n_qubits, min_abs_coeff = 1e-10, num_shots = 1, seed = None))]
            pub fn sample(
                py: Python<'_>,
                prog: &crate::stim_program::PyStimProgram,
                n_qubits: usize,
                min_abs_coeff: f64,
                num_shots: usize,
                seed: Option<u64>,
            ) -> pyo3::PyResult<Vec<Vec<u8>>> {
                // `prog` was already validated at `StimProgram.parse()` time;
                // use the validated path to skip redundant re-validation.
                let raw = py.detach(|| {
                    ppvm_stim::sample_validated(
                        &prog.0.instructions,
                        prog.0.measurement_count(),
                        num_shots,
                        |i| match seed {
                            Some(s) => <$type>::new_with_seed(
                                n_qubits,
                                min_abs_coeff,
                                s.wrapping_add(i as u64),
                            ),
                            None => <$type>::new(n_qubits, min_abs_coeff),
                        },
                    )
                });
                Ok(raw
                    .into_iter()
                    .map(|shot| {
                        shot.into_iter()
                            .map(crate::interface_tableau::measurement_to_u8)
                            .collect()
                    })
                    .collect())
            }

            /// Fork this tableau, cloning all quantum state but reinitializing the RNG.
            /// If `seed` is provided, the new RNG is seeded deterministically; otherwise
            /// it is seeded from OS entropy, giving an independent random sequence.
            ///
            /// Use this when branching a simulation into independent trajectories.
            /// To preserve the RNG state exactly (e.g. for checkpointing), use
            /// `copy.copy()` or `copy.deepcopy()` instead.
            #[pyo3(signature = (seed = None))]
            pub fn fork(&self, seed: Option<u64>) -> Self {
                Self {
                    inner: self.inner.fork(seed),
                }
            }

            /// Return a shallow copy of this tableau, including its RNG state.
            ///
            /// Both the original and the copy will produce identical random sequences
            /// from this point forward. To get an independent copy with a fresh RNG,
            /// use `fork()` instead.
            fn __copy__(&self) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }

            /// Return a deep copy of this tableau, including its RNG state.
            ///
            /// Both the original and the copy will produce identical random sequences
            /// from this point forward. To get an independent copy with a fresh RNG,
            /// use `fork()` instead.
            fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }
        }
    };
}
