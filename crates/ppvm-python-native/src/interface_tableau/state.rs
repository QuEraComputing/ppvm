// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

#[cfg(feature = "legacy")]
macro_rules! expectation_word {
    ($storage:expr, $word:expr) => {{
        let word: PauliWord<[usize; $storage]> = $word.into();
        word
    }};
}

#[cfg(feature = "traits-2")]
macro_rules! expectation_word {
    ($storage:expr, $word:expr) => {{
        let word: ppvm_pauli_sum_2::PauliWord<[usize; $storage]> = $word.into();
        word
    }};
}
/// Construct the wrapped tableau. Legacy seeds the tableau's embedded RNG.
#[cfg(feature = "legacy")]
macro_rules! new_tableau {
    ($type:ty, $n:expr, $min_abs_coeff:expr, $seed:expr) => {{
        let tab: $type = match $seed {
            Some(s) => <$type>::new_with_seed($n, $min_abs_coeff, s),
            None => <$type>::new($n, $min_abs_coeff),
        };
        Self { inner: tab }
    }};
}

/// Construct the wrapped tableau. Under `-2` the tableau is pure state and the
/// seed goes to the wrapper's own RNG.
#[cfg(feature = "traits-2")]
macro_rules! new_tableau {
    ($type:ty, $n:expr, $min_abs_coeff:expr, $seed:expr) => {
        wrap!(<$type>::new($n, $min_abs_coeff), $seed)
    };
}

macro_rules! create_tableau_state {
    ($name: ident, $type: ident, $storage:expr, $indexType: ident) => {
        #[pyclass]
        pub struct $name {
            inner: $type,
            /// See `§ Where the randomness lives` in `crate::backend`.
            #[cfg(feature = "traits-2")]
            rng: rand::rngs::SmallRng,
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, seed = None))]
            pub fn new(n_qubits: usize, min_abs_coeff: f64, seed: Option<u64>) -> Self {
                new_tableau!($type, n_qubits, min_abs_coeff, seed)
            }

            fn __repr__(&self) -> String {
                // TODO: expose some more details e.g. for debugging
                format!("{}", self.inner)
            }

            fn __str__(&self) -> String {
                self.inner.to_string()
            }

            pub fn measure(&mut self, addr0: usize) -> i64 {
                measurement_to_u8(draw!(self.inner.measure(addr0))) as i64
            }

            pub fn measure_many(&mut self, targets: Vec<usize>) -> Vec<i64> {
                draw!(self.inner.measure_many(targets.as_slice()))
                    .into_iter()
                    .map(|m| measurement_to_u8(m) as i64)
                    .collect()
            }

            pub fn current_measurement_record(&self) -> Vec<i64> {
                self.inner
                    .current_measurement_record()
                    .iter()
                    .map(|m| measurement_to_u8(*m) as i64)
                    .collect()
            }

            /// Snapshot of the sparse coefficient vector as `{index: amplitude}`.
            ///
            /// Keys are basis-state indices (Python ints, lossless at every
            /// width); values are complex amplitudes. This is a copy — mutating
            /// it does not touch the tableau's internal state.
            pub fn coefficients<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                let dict = PyDict::new(py);
                // `int(str)` handles index widths beyond u128 (bnum types).
                let int_ctor = py.import("builtins")?.getattr("int")?;
                for (coeff, idx) in self.inner.coefficients.iter() {
                    let key = int_ctor.call1((idx.to_string(),))?;
                    let value = PyComplex::from_doubles(py, coeff.re, coeff.im);
                    dict.set_item(key, value)?;
                }
                Ok(dict)
            }

            /// Number of branches stored in the coefficient vector.
            pub fn num_coefficients(&self) -> usize {
                self.inner.coefficients.len()
            }

            /// `⟨ψ|word|ψ⟩` for the multi-qubit Pauli string `word`.
            pub fn expectation(&self, word: String) -> f64 {
                let w = expectation_word!($storage, word);
                self.inner.expectation(&w)
            }

            /// `Σ_{P matches pattern} ⟨ψ|P|ψ⟩`.
            pub fn trace(&self, pattern: String) -> f64 {
                let pat: PauliPattern = pattern.into();
                self.inner.trace(&pat)
            }
        }
    };
}
