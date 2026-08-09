// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// Run a validated program against the wrapped tableau.
#[cfg(feature = "legacy")]
macro_rules! execute_program {
    ($this:ident, $instructions:expr, $results:expr) => {
        ppvm_stim::execute_validated($instructions, &mut $this.inner, $results)
    };
}

/// Run a validated program against the wrapped tableau, feeding it the
/// wrapper's RNG.
#[cfg(feature = "traits-2")]
macro_rules! execute_program {
    ($this:ident, $instructions:expr, $results:expr) => {
        ppvm_stim::execute_validated_with_rng(
            $instructions,
            &mut $this.inner,
            $results,
            &mut $this.rng,
        )
    };
}

/// Build the per-shot state for `sample`.
///
/// Legacy folds the shot seed into each fresh tableau; `-2` hands it to the
/// per-shot RNG instead. Either way shot `i` derives from `seed + i`, so the
/// two backends produce the same shots for the same seed.
#[cfg(feature = "legacy")]
macro_rules! sample_shots {
    ($type:ty, $instructions:expr, $count:expr, $num_shots:expr, $n:expr, $min:expr, $seed:expr) => {
        ppvm_stim::sample_validated($instructions, $count, $num_shots, |i| match $seed {
            Some(s) => <$type>::new_with_seed($n, $min, s.wrapping_add(i as u64)),
            None => <$type>::new($n, $min),
        })
    };
}

/// Build the per-shot state and RNG for `sample`.
#[cfg(feature = "traits-2")]
macro_rules! sample_shots {
    ($type:ty, $instructions:expr, $count:expr, $num_shots:expr, $n:expr, $min:expr, $seed:expr) => {
        ppvm_stim::sample_validated_with_rng(
            $instructions,
            $count,
            $num_shots,
            |_| <$type>::new($n, $min),
            |i| crate::backend::make_rng($seed.map(|s: u64| s.wrapping_add(i as u64))),
        )
    };
}

/// Clone the state and reseed. Legacy's `fork` reseeds the tableau itself.
#[cfg(feature = "legacy")]
macro_rules! fork_tableau {
    ($this:ident, $seed:expr) => {
        Self {
            inner: $this.inner.fork($seed),
        }
    };
}

/// Clone the state and reseed. Under `-2` the tableau is pure state, so the
/// reseed is entirely the wrapper's business.
#[cfg(feature = "traits-2")]
macro_rules! fork_tableau {
    ($this:ident, $seed:expr) => {
        wrap!($this.inner.fork(), $seed)
    };
}

macro_rules! create_tableau_stim {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pymethods]
        impl $name {
            pub fn run(
                &mut self,
                prog: &crate::stim_program::PyStimProgram,
            ) -> pyo3::PyResult<Vec<u8>> {
                let mut results = Vec::with_capacity(prog.measurement_count());
                execute_program!(self, &prog.instructions, &mut results);
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
                    sample_shots!(
                        $type,
                        &prog.0.instructions,
                        prog.0.measurement_count(),
                        num_shots,
                        n_qubits,
                        min_abs_coeff,
                        seed
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
                fork_tableau!(self, seed)
            }

            /// Return a shallow copy of this tableau, including its RNG state.
            ///
            /// Both the original and the copy will produce identical random sequences
            /// from this point forward. To get an independent copy with a fresh RNG,
            /// use `fork()` instead.
            fn __copy__(&self) -> Self {
                wrap_cloned!(self)
            }

            /// Return a deep copy of this tableau, including its RNG state.
            ///
            /// Both the original and the copy will produce identical random sequences
            /// from this point forward. To get an independent copy with a fresh RNG,
            /// use `fork()` instead.
            fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
                wrap_cloned!(self)
            }
        }
    };
}
