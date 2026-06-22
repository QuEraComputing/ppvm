// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bnum::types::{U256, U512, U1024, U2048};
use paste::paste;
use ppvm_tableau::prelude::*;
use pyo3::prelude::*;

pub(crate) fn measurement_to_u8(m: Option<bool>) -> u8 {
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

            pub fn measure(&mut self, addr0: usize) -> i64 {
                measurement_to_u8(self.inner.measure(addr0)) as i64
            }

            pub fn measure_many(&mut self, targets: Vec<usize>) -> Vec<i64> {
                self.inner
                    .measure_many(targets.as_slice())
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

            // clifford
            pub fn x(&mut self, targets: Vec<usize>) {
                self.inner.x(targets.as_slice());
            }

            pub fn y(&mut self, targets: Vec<usize>) {
                self.inner.y(targets.as_slice());
            }

            pub fn z(&mut self, targets: Vec<usize>) {
                self.inner.z(targets.as_slice());
            }

            pub fn h(&mut self, targets: Vec<usize>) {
                self.inner.h(targets.as_slice());
            }

            pub fn s(&mut self, targets: Vec<usize>) {
                self.inner.s(targets.as_slice());
            }

            pub fn s_dag(&mut self, targets: Vec<usize>) {
                self.inner.s_dag(targets.as_slice());
            }

            // clifford extensions
            pub fn sqrt_x(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_x(targets.as_slice());
            }

            pub fn sqrt_x_dag(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_x_dag(targets.as_slice());
            }

            pub fn sqrt_y(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_y(targets.as_slice());
            }

            pub fn sqrt_y_dag(&mut self, targets: Vec<usize>) {
                self.inner.sqrt_y_dag(targets.as_slice());
            }

            pub fn t(&mut self, targets: Vec<usize>) {
                self.inner.t(targets.as_slice());
            }

            pub fn t_dag(&mut self, targets: Vec<usize>) {
                self.inner.t_dag(targets.as_slice());
            }

            // two-qubit clifford (+ stim aliases)
            pub fn cnot(&mut self, targets: Vec<usize>) {
                self.inner.cnot(targets.as_slice());
            }

            pub fn cx(&mut self, targets: Vec<usize>) {
                self.inner.cnot(targets.as_slice());
            }

            pub fn zcx(&mut self, targets: Vec<usize>) {
                self.inner.cnot(targets.as_slice());
            }

            pub fn cy(&mut self, targets: Vec<usize>) {
                self.inner.cy(targets.as_slice());
            }

            pub fn zcy(&mut self, targets: Vec<usize>) {
                self.inner.cy(targets.as_slice());
            }

            pub fn cz(&mut self, targets: Vec<usize>) {
                self.inner.cz(targets.as_slice());
            }

            pub fn zcz(&mut self, targets: Vec<usize>) {
                self.inner.cz(targets.as_slice());
            }

            // rot1
            pub fn rx(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rx(targets.as_slice(), theta);
            }

            pub fn ry(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.ry(targets.as_slice(), theta);
            }

            pub fn rz(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rz(targets.as_slice(), theta);
            }

            pub fn u3(&mut self, addr0: usize, theta: f64, phi: f64, lam: f64) {
                self.inner.u3(addr0, theta, phi, lam);
            }

            // rot2
            pub fn rxx(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rxx(targets.as_slice(), theta);
            }

            pub fn ryy(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.ryy(targets.as_slice(), theta);
            }

            pub fn rzz(&mut self, targets: Vec<usize>, theta: f64) {
                self.inner.rzz(targets.as_slice(), theta);
            }

            // noise
            pub fn x_error(&mut self, targets: Vec<usize>, p: f64) {
                self.inner.x_error(targets.as_slice(), p);
            }

            pub fn y_error(&mut self, targets: Vec<usize>, p: f64) {
                self.inner.y_error(targets.as_slice(), p);
            }

            pub fn z_error(&mut self, targets: Vec<usize>, p: f64) {
                self.inner.z_error(targets.as_slice(), p);
            }

            pub fn pauli_error(&mut self, targets: Vec<usize>, p: [f64; 3]) {
                self.inner.pauli_error(targets.as_slice(), p);
            }

            pub fn depolarize1(&mut self, targets: Vec<usize>, p: f64) {
                self.inner.depolarize1(targets.as_slice(), p);
            }

            pub fn depolarize2(&mut self, targets: Vec<usize>, p: f64) {
                self.inner.depolarize2(targets.as_slice(), p);
            }

            pub fn two_qubit_pauli_error(&mut self, targets: Vec<usize>, p: [f64; 15]) {
                self.inner.two_qubit_pauli_error(targets.as_slice(), p);
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

            pub fn reset(&mut self, targets: Vec<usize>) {
                self.inner.reset(targets.as_slice());
            }

            pub fn reset_x(&mut self, targets: Vec<usize>) {
                self.inner.reset_x(targets.as_slice());
            }

            pub fn reset_y(&mut self, targets: Vec<usize>) {
                self.inner.reset_y(targets.as_slice());
            }

            pub fn reset_z(&mut self, targets: Vec<usize>) {
                self.inner.reset_z(targets.as_slice());
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
                            Some(s) => GeneralizedTableau::<$type, $indexType>::new_with_seed(
                                n_qubits,
                                min_abs_coeff,
                                s.wrapping_add(i as u64),
                            ),
                            None => GeneralizedTableau::<$type, $indexType>::new(
                                n_qubits,
                                min_abs_coeff,
                            ),
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
