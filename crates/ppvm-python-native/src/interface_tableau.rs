// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use paste::paste;
use ppvm_tableau::prelude::*;
use pyo3::prelude::*;

fn measurement_to_u8(m: Option<bool>) -> u8 {
    match m {
        Some(false) => 0,
        Some(true) => 1,
        None => 2,
    }
}

macro_rules! create_interface {
    ($name: ident, $type: ident, $indexType: ident) => {
        #[pyclass]
        pub struct $name {
            inner: GeneralizedTableau<$type, $indexType>,
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, seed = None))]
            pub fn new(n_qubits: usize, min_abs_coeff: f64, seed: Option<u64>) -> Self {
                let tab: GeneralizedTableau<$type, $indexType> = match seed {
                    Some(s) => GeneralizedTableau::new_with_seed(n_qubits, min_abs_coeff, s),
                    None => GeneralizedTableau::new(n_qubits, min_abs_coeff),
                };
                Self { inner: tab }
            }

            fn __repr__(&self) -> String {
                // TODO: expose some more details e.g. for debugging
                format!("{}", self.inner)
            }

            fn __str__(&self) -> String {
                self.inner.to_string()
            }

            pub fn measure(&mut self, addr0: usize) -> Option<bool> {
                self.inner.measure(addr0)
            }

            // clifford
            pub fn x(&mut self, addr0: usize) {
                self.inner.x(addr0);
            }

            pub fn y(&mut self, addr0: usize) {
                self.inner.y(addr0);
            }

            pub fn z(&mut self, addr0: usize) {
                self.inner.z(addr0);
            }

            pub fn h(&mut self, addr0: usize) {
                self.inner.h(addr0);
            }

            pub fn s(&mut self, addr0: usize) {
                self.inner.s(addr0);
            }

            pub fn s_adj(&mut self, addr0: usize) {
                self.inner.s_adj(addr0);
            }

            // clifford extensions
            pub fn sqrt_x(&mut self, addr0: usize) {
                self.inner.sqrt_x(addr0);
            }

            pub fn sqrt_x_adj(&mut self, addr0: usize) {
                self.inner.sqrt_x_adj(addr0);
            }

            pub fn sqrt_y(&mut self, addr0: usize) {
                self.inner.sqrt_y(addr0);
            }

            pub fn sqrt_y_adj(&mut self, addr0: usize) {
                self.inner.sqrt_y_adj(addr0);
            }

            pub fn cnot(&mut self, addr0: usize, addr1: usize) {
                self.inner.cnot(addr0, addr1);
            }

            pub fn cy(&mut self, addr0: usize, addr1: usize) {
                self.inner.cy(addr0, addr1);
            }

            pub fn cz(&mut self, addr0: usize, addr1: usize) {
                self.inner.cz(addr0, addr1);
            }

            pub fn t(&mut self, addr0: usize) {
                self.inner.t(addr0);
            }

            pub fn t_adj(&mut self, addr0: usize) {
                self.inner.t_adj(addr0);
            }

            // rot1
            pub fn rx(&mut self, addr0: usize, theta: f64) {
                self.inner.rx(addr0, theta);
            }

            pub fn ry(&mut self, addr0: usize, theta: f64) {
                self.inner.ry(addr0, theta);
            }

            pub fn rz(&mut self, addr0: usize, theta: f64) {
                self.inner.rz(addr0, theta);
            }

            pub fn u3(&mut self, addr0: usize, theta: f64, phi: f64, lam: f64) {
                self.inner.u3(addr0, theta, phi, lam);
            }

            // rot2
            pub fn rxx(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rxx(addr0, addr1, theta);
            }

            pub fn ryy(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.ryy(addr0, addr1, theta);
            }

            pub fn rzz(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rzz(addr0, addr1, theta);
            }

            // noise
            pub fn pauli_error(&mut self, addr0: usize, p: [f64; 3]) {
                self.inner.pauli_error(addr0, p);
            }

            pub fn depolarize(&mut self, addr0: usize, p: f64) {
                self.inner.depolarize(addr0, p);
            }

            pub fn depolarize2(&mut self, addr0: usize, addr1: usize, p: f64) {
                self.inner.depolarize2(addr0, addr1, p);
            }

            pub fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [f64; 15]) {
                self.inner.two_qubit_pauli_error(addr0, addr1, p);
            }

            pub fn loss_channel(&mut self, addr0: usize, p: f64) {
                self.inner.loss_channel(addr0, p);
            }

            pub fn correlated_loss_channel(&mut self, addr0: usize, addr1: usize, p: [f64; 3]) {
                self.inner.correlated_loss_channel(addr0, addr1, p);
            }

            pub fn reset_loss_channel(&mut self, addr0: usize) {
                self.inner.reset_loss_channel(addr0);
            }

            pub fn reset(&mut self, addr0: usize) {
                self.inner.reset(addr0);
            }

            pub fn is_lost(&self, addr0: usize) -> bool {
                self.inner.is_lost[addr0]
            }

            pub fn loss_values(&self) -> Vec<bool> {
                self.inner.is_lost.clone()
            }

            pub fn run(
                &mut self,
                prog: &crate::stim_program::PyStimProgram,
            ) -> pyo3::PyResult<Vec<u8>> {
                let mut results = Vec::with_capacity(prog.measurement_count());
                ppvm_stim::execute_prepared(&prog.instructions, &mut self.inner, &mut results);
                Ok(results
                    .into_iter()
                    .map(crate::interface_tableau::measurement_to_u8)
                    .collect())
            }

            /// Multi-shot sampling: builds a fresh tableau per shot.
            #[staticmethod]
            #[pyo3(signature = (prog, n_qubits, min_abs_coeff = 1e-10, num_shots = 1, seed = None))]
            pub fn sample(
                prog: &crate::stim_program::PyStimProgram,
                n_qubits: usize,
                min_abs_coeff: f64,
                num_shots: usize,
                seed: Option<u64>,
            ) -> pyo3::PyResult<Vec<Vec<u8>>> {
                let mut next_seed = seed;
                let count = prog.measurement_count();
                Ok((0..num_shots)
                    .map(|_| {
                        let s = next_seed;
                        if let Some(ref mut v) = next_seed {
                            *v = v.wrapping_add(1);
                        }
                        let mut tab = match s {
                            Some(s) => GeneralizedTableau::<$type, $indexType>::new_with_seed(
                                n_qubits,
                                min_abs_coeff,
                                s,
                            ),
                            None => GeneralizedTableau::<$type, $indexType>::new(
                                n_qubits,
                                min_abs_coeff,
                            ),
                        };
                        let mut shot = Vec::with_capacity(count);
                        ppvm_stim::execute_prepared(&prog.instructions, &mut tab, &mut shot);
                        shot
                    })
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

macro_rules! create_interface_range {
    ($name: ident, $indexType: ident, $( $n: expr),* ) => {
        paste! {
        $(
            type [<$name$n>] = config::fx64hash::Byte8F64<$n>;
            // type [<$name$n>] = config::indexmap::ByteFxHashF64<$n>;
            create_interface!([<GeneralizedTableau$n>], [<$name$n>], $indexType);
        )*
    }
    };
}

// up to 64 qubits
create_interface_range!(IndexMapFxHash, usize, 1);

// 64 - 128 qubits
create_interface_range!(IndexMapFxHash, u128, 2);

// 128 - 256 qubits
create_interface_range!(IndexMapFxHash, U256, 3, 4);

create_interface_range!(IndexMapFxHash, U512, 5, 6, 7, 8);

create_interface_range!(IndexMapFxHash, U1024, 9, 10, 11, 12, 13, 14, 15, 16);

create_interface_range!(
    IndexMapFxHash,
    U2048,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32
);
