// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use paste::paste;
use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::{
    CoefficientThreshold, CombinedStrategy, MaxLossWeight, MaxPauliWeight,
};
use pyo3::prelude::*;

macro_rules! create_interface_loss_methods {
    ($name: ident, $type: ident, false) => {};
    ($name: ident, $type: ident, true) => {
        #[pymethods]
        impl $name {
            pub fn loss_channel(&mut self, addr0: usize, p: f64) {
                self.inner.loss_channel(addr0, p);
                self.inner.truncate();
            }

            pub fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [f64; 3]) {
                self.inner.correlated_loss_channel(addr0, addr1, p);
                self.inner.truncate();
            }

            pub fn reset_loss_channel(&mut self, addr0: usize) {
                self.inner.reset_loss_channel(addr0);
                self.inner.truncate();
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
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, max_pauli_weight = usize::MAX, max_loss_weight = usize::MAX, terms = Vec::<String>::new(), coefficients = Vec::<f64>::new()))]
            pub fn new(
                n_qubits: usize,
                min_abs_coeff: f64,
                max_pauli_weight: usize,
                max_loss_weight: usize,
                terms: Vec<String>,
                coefficients: Vec<f64>
            ) -> Self {
                let _ = max_loss_weight; // unused in non-loss variants
                let strategy = create_strategy!($loss, min_abs_coeff, max_pauli_weight, max_loss_weight);
                let mut ps = PauliSum::<$type>::builder()
                    .n_qubits(n_qubits)
                    .strategy(strategy)
                    .capacity(n_qubits)
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

            // clifford
            pub fn x(&mut self, addr0: usize) {
                self.inner.x(addr0);
                self.inner.truncate();
            }

            pub fn y(&mut self, addr0: usize) {
                self.inner.y(addr0);
                self.inner.truncate();
            }

            pub fn z(&mut self, addr0: usize) {
                self.inner.z(addr0);
                self.inner.truncate();
            }

            pub fn h(&mut self, addr0: usize) {
                self.inner.h(addr0);
                self.inner.truncate();
            }

            pub fn s(&mut self, addr0: usize) {
                self.inner.s(addr0);
                self.inner.truncate();
            }

            pub fn cnot(&mut self, addr0: usize, addr1: usize) {
                self.inner.cnot(addr0, addr1);
                self.inner.truncate();
            }

            pub fn cz(&mut self, addr0: usize, addr1: usize) {
                self.inner.cz(addr0, addr1);
                self.inner.truncate();
            }

            // clifford extensions
            pub fn s_adj(&mut self, addr0: usize) {
                self.inner.s_adj(addr0);
                self.inner.truncate();
            }

            pub fn sqrt_x(&mut self, addr0: usize) {
                self.inner.sqrt_x(addr0);
                self.inner.truncate();
            }

            pub fn sqrt_y(&mut self, addr0: usize) {
                self.inner.sqrt_y(addr0);
                self.inner.truncate();
            }

            pub fn sqrt_x_adj(&mut self, addr0: usize) {
                self.inner.sqrt_x_adj(addr0);
                self.inner.truncate();
            }

            pub fn sqrt_y_adj(&mut self, addr0: usize) {
                self.inner.sqrt_y_adj(addr0);
                self.inner.truncate();
            }

            // rot1
            pub fn rx(&mut self, addr0: usize, theta: f64) {
                self.inner.rx(addr0, theta);
                self.inner.truncate();
            }

            pub fn ry(&mut self, addr0: usize, theta: f64) {
                self.inner.ry(addr0, theta);
                self.inner.truncate();
            }

            pub fn rz(&mut self, addr0: usize, theta: f64) {
                self.inner.rz(addr0, theta);
                self.inner.truncate();
            }

            // rot2
            pub fn rxx(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rxx(addr0, addr1, theta);
                self.inner.truncate();
            }

            pub fn ryy(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.ryy(addr0, addr1, theta);
                self.inner.truncate();
            }

            pub fn rzz(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rzz(addr0, addr1, theta);
                self.inner.truncate();
            }

            // U(1)-conserving exchange/Heisenberg-style gates
            pub fn exchange(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.exchange(addr0, addr1, theta);
                self.inner.truncate();
            }

            pub fn xyzz(&mut self, addr0: usize, addr1: usize, theta_xy: f64, theta_zz: f64) {
                self.inner.xyzz(addr0, addr1, theta_xy, theta_zz);
                self.inner.truncate();
            }

            // noise
            pub fn pauli_error(&mut self, addr0: usize, p: [f64; 3]) {
                self.inner.pauli_error(addr0, p);
                self.inner.truncate();
            }

            pub fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [f64; 15]) {
                self.inner.two_qubit_pauli_error(addr0, addr1, p);
                self.inner.truncate();
            }

            pub fn depolarize(&mut self, addr0: usize, p: f64) {
                self.inner.depolarize(addr0, p);
                self.inner.truncate();
            }

            pub fn depolarize2(&mut self, addr0: usize, addr1: usize, p: f64) {
                self.inner.depolarize2(addr0, addr1, p);
                self.inner.truncate();
            }

            pub fn amplitude_damping(&mut self, addr0: usize, gamma: f64) {
                self.inner.amplitude_damping(addr0, gamma);
                self.inner.truncate();
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
