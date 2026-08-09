// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! create_interface_state {
    ($name: ident, $type: ident, $loss: tt) => {
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
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, max_pauli_weight = usize::MAX, max_loss_weight = usize::MAX, terms = Vec::<String>::new(), coefficients = Vec::<f64>::new(), preserve_strings = Vec::<String>::new()))]
            // The atomic digest cache is excluded from structural Eq/Hash and is
            // invalidated only through `&mut` word operations. Stored keys cannot
            // change map identity; Clippy cannot see that invariant.
            #[allow(clippy::mutable_key_type)]
            pub fn new(
                n_qubits: usize,
                min_abs_coeff: f64,
                max_pauli_weight: usize,
                max_loss_weight: usize,
                terms: Vec<String>,
                coefficients: Vec<f64>,
                preserve_strings: Vec<String>,
            ) -> Self {
                let _ = max_loss_weight; // unused in non-loss variants
                let strategy = create_strategy!($loss, min_abs_coeff, max_pauli_weight, max_loss_weight);
                assert_eq!(
                    terms.len(),
                    coefficients.len(),
                    "Initial terms and coefficients need to be of same length!"
                );

                let ps = construct_pauli_sum!(
                    $type,
                    strategy,
                    n_qubits,
                    preserve_strings,
                    terms,
                    coefficients
                );

                wrap!(ps, None)
            }

            fn __repr__(&self) -> String {
                // TODO: expose some more details e.g. for debugging
                format!("PauliSum({})", self.inner)
            }

            fn __str__(&self) -> String {
                self.inner.to_string()
            }

            pub fn trace(&self, pattern: String) -> f64 {
                let pat: PauliPattern = pattern.into();
                let result = self.inner.trace(&pat);
                result
            }

            pub fn overlap_with_zero(&self) -> f64 {
                self.trace("Z?*".to_owned())
            }

            pub fn overlap(&self, other: &Self) -> f64 {
                self.inner.overlap(&other.inner)
            }

            // NOTE: macros can't be used in pymethods block
            // could either use multiple-pymethods feature (adds dependencies)
            // or better yet create working impl for all strategies
            //
            // Every gate exposes a `truncate: bool = True` kwarg. When
            // `True` (the default) the inner `truncate()` runs immediately
            // after the gate, matching the historical behaviour. Pass
            // `truncate=False` to defer the cut — useful for chaining
            // commuting gates (e.g. `rxx + ryy` on the same edge, or any
            // other U(1)/Z₂-conserving composition) where truncating
            // between them would break a conserved-charge structure that
            // truncating only once at the end preserves.

            /// Explicit truncate. Use with `truncate=False` on the gates
            /// above to control exactly when the active strategy fires.
            pub fn truncate(&mut self) {
                self.inner.truncate();
            }
        }
    };
}
