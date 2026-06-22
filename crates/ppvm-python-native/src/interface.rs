// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashSet;

use paste::paste;
use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold, CombinedStrategy, MaxLossWeight, MaxPauliWeight,
};
use pyo3::prelude::*;

macro_rules! create_interface_loss_methods {
    ($name: ident, $type: ident, false) => {};
    ($name: ident, $type: ident, true) => {
        #[pymethods]
        impl $name {
            #[pyo3(signature = (addr0, p, truncate = true))]
            pub fn loss_channel(&mut self, addr0: usize, p: f64, truncate: bool) {
                self.inner.loss_channel(addr0, p);
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (addr0, addr1, p, truncate = true))]
            pub fn correlated_loss_channel(
                &mut self,
                addr0: usize,
                addr1: usize,
                p: [f64; 3],
                truncate: bool,
            ) {
                self.inner.correlated_loss_channel(addr0, addr1, p);
                if truncate {
                    self.inner.truncate();
                }
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn reset_loss_channel(&mut self, addr0: usize, truncate: bool) {
                self.inner.reset_loss_channel(addr0);
                if truncate {
                    self.inner.truncate();
                }
            }
        }
    };
}

macro_rules! create_strategy {
    (false, $min_abs_coeff:ident, $max_pauli_weight:ident, $_max_loss_weight:ident) => {
        CombinedStrategy(
            CoefficientThreshold($min_abs_coeff),
            MaxPauliWeight($max_pauli_weight),
        )
    };
    (true, $min_abs_coeff:ident, $max_pauli_weight:ident, $max_loss_weight:ident) => {
        CombinedStrategy(
            CombinedStrategy(
                CoefficientThreshold($min_abs_coeff),
                MaxPauliWeight($max_pauli_weight),
            ),
            MaxLossWeight($max_loss_weight),
        )
    };
}

// adapted from https://pyo3.rs/v0.27.1/class.html#no-generic-parameters
macro_rules! create_interface {
    ($name: ident, $type: ident, $loss: tt) => {
        #[pyclass]
        pub struct $name {
            inner: PauliSum<$type>,
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, max_pauli_weight = usize::MAX, max_loss_weight = usize::MAX, terms = Vec::<String>::new(), coefficients = Vec::<f64>::new(), preserve_strings = Vec::<String>::new()))]
            #[allow(clippy::too_many_arguments)]
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
                let preserve_set: HashSet<_> =
                    preserve_strings.into_iter().map(Into::into).collect();
                let mut ps = PauliSum::<$type>::builder()
                    .n_qubits(n_qubits)
                    .strategy(strategy)
                    .capacity(n_qubits)
                    .preserve_strings(preserve_set)
                    .build();

                assert_eq!(
                    terms.len(),
                    coefficients.len(),
                    "Initial terms and coefficients need to be of same length!"
                );

                for (term, c) in terms.iter().zip(coefficients.iter()) {
                    ps += (term.to_owned(), c.to_owned());
                }

                Self { inner: ps }
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

            // clifford
            #[pyo3(signature = (targets, truncate = true))]
            pub fn x(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.x_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn y(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.y_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn z(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.z_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn h(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.h_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn s(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.s_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            // two-qubit clifford (+ stim aliases)
            #[pyo3(signature = (targets, truncate = true))]
            pub fn cnot(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cnot_batch(&pairs);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cx(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cnot(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcx(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cnot(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cz(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cz_batch(&pairs);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcz(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cz(targets, truncate)
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn cy(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.cy_batch(&pairs);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn zcy(&mut self, targets: Vec<usize>, truncate: bool) -> PyResult<()> {
                self.cy(targets, truncate)
            }

            // clifford extensions
            #[pyo3(signature = (targets, truncate = true))]
            pub fn s_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.s_adj_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_x(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_x_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_y(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_y_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_x_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_x_adj_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, truncate = true))]
            pub fn sqrt_y_dag(&mut self, targets: Vec<usize>, truncate: bool) {
                self.inner.sqrt_y_adj_batch(targets.as_slice());
                if truncate { self.inner.truncate(); }
            }

            // rot1
            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rx(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.rx_batch(targets.as_slice(), theta);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn ry(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.ry_batch(targets.as_slice(), theta);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rz(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) {
                self.inner.rz_batch(targets.as_slice(), theta);
                if truncate { self.inner.truncate(); }
            }

            // rot2
            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rxx(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rxx_batch(&pairs, theta);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn ryy(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.ryy_batch(&pairs, theta);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, theta, truncate = true))]
            pub fn rzz(&mut self, targets: Vec<usize>, theta: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.rzz_batch(&pairs, theta);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            // noise
            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn x_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                self.inner.x_error_batch(targets.as_slice(), p);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn y_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                self.inner.y_error_batch(targets.as_slice(), p);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn z_error(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                self.inner.z_error_batch(targets.as_slice(), p);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn pauli_error(&mut self, targets: Vec<usize>, p: [f64; 3], truncate: bool) {
                self.inner.pauli_error_batch(targets.as_slice(), p);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn two_qubit_pauli_error(
                &mut self,
                targets: Vec<usize>,
                p: [f64; 15],
                truncate: bool,
            ) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.two_qubit_pauli_error_batch(&pairs, p);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn depolarize1(&mut self, targets: Vec<usize>, p: f64, truncate: bool) {
                self.inner.depolarize1_batch(targets.as_slice(), p);
                if truncate { self.inner.truncate(); }
            }

            #[pyo3(signature = (targets, p, truncate = true))]
            pub fn depolarize2(&mut self, targets: Vec<usize>, p: f64, truncate: bool) -> PyResult<()> {
                let pairs = crate::flat_pairs(&targets)?;
                self.inner.depolarize2_batch(&pairs, p);
                if truncate { self.inner.truncate(); }
                Ok(())
            }

            #[pyo3(signature = (addr0, gamma, truncate = true))]
            pub fn amplitude_damping(&mut self, addr0: usize, gamma: f64, truncate: bool) {
                self.inner.amplitude_damping(addr0, gamma);
                if truncate { self.inner.truncate(); }
            }


            // some python niceties

            fn __copy__(&self) -> Self {
                Self { inner: self.inner.clone() }
            }

            fn __richcmp__(&self, other: PyRef<$name>, op: pyo3::basic::CompareOp) -> PyResult<bool> {
                match op {
                    pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
                    pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
                    _ => Err(pyo3::exceptions::PyNotImplementedError::new_err("Only equality and inequality comparisons are supported for PauliSum.")),
                }
            }

            fn __deepcopy__(&self, _memo: &Bound<'_, PyAny>) -> Self {
                Self { inner: self.inner.clone() }
            }

            fn __len__(&self) -> usize {
                self.inner.len()
            }

            pub fn terms(&self) -> Vec<(String, f64)> {
                self.inner.data().iter().map(|(k, v)| (k.to_string(), *v)).collect()
            }

            pub fn weights(&self) -> Vec<(String, usize)> {
                self.inner.data().iter().map(|(k, _v)| (k.to_string(), k.weight())).collect()
            }

            pub fn current_max_weight(&self) -> usize {
                self.inner.data().iter().map(|(k, _v)| k.weight()).max().unwrap_or(0)
            }
        }

        create_interface_loss_methods!($name, $type, $loss);
    };
}

macro_rules! create_interface_range {
    ($name: ident, false, $( $n: expr),* ) => {
        paste! {
        $(
            type [<$name$n>] = config::indexmap::ByteFxHashF64<{(2 as usize).pow($n)}, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
            create_interface!([<PauliSum$name$n>], [<$name$n>], false);
        )*
    }
    };

    ($name: ident, true, $( $n: expr),* ) => {
        paste! {
        $(
            type [<Loss$name$n>] = config::indexmap::ByteFxHashF64<{(2 as usize).pow($n)}, CombinedStrategy<CombinedStrategy<CoefficientThreshold, MaxPauliWeight>, MaxLossWeight>, LossyPauliWord<[u8; {(2 as usize).pow($n)}]>>;
            create_interface!([<PauliSumLoss$name$n>], [<Loss$name$n>], true);
        )*
    }
    };
}

create_interface_range!(
    IndexMapFxHash,
    false,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15
);

create_interface_range!(
    IndexMapFxHash,
    true,
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15
);
