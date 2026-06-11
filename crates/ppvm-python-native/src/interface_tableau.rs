// SPDX-FileCopyrightText: 2026 The PPVM Authors
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

            // All gate methods take a `truncate: bool = true` kwarg
            // purely for API symmetry with `PauliSum` (which uses it to
            // defer the auto-truncate). `GeneralizedTableau` has no
            // analogous truncation step — the tableau representation is
            // exact — so the kwarg is silently ignored here. Keeping the
            // signature parallel lets the Python `RotationsMixin` /
            // `CliffordMixin` / etc. work uniformly across both backends.

            // clifford
            #[pyo3(signature = (addr0, truncate = true))]
            pub fn x(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.x(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn y(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.y(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn z(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.z(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn h(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.h(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn s(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.s(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn s_adj(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.s_adj(addr0);
            }

            // clifford extensions
            #[pyo3(signature = (addr0, truncate = true))]
            pub fn sqrt_x(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.sqrt_x(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn sqrt_x_adj(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.sqrt_x_adj(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn sqrt_y(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.sqrt_y(addr0);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn sqrt_y_adj(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
                self.inner.sqrt_y_adj(addr0);
            }

            #[pyo3(signature = (addr0, addr1, truncate = true))]
            pub fn cnot(&mut self, addr0: usize, addr1: usize, truncate: bool) {
                let _ = truncate;
                self.inner.cnot(addr0, addr1);
            }

            #[pyo3(signature = (addr0, addr1, truncate = true))]
            pub fn cy(&mut self, addr0: usize, addr1: usize, truncate: bool) {
                let _ = truncate;
                self.inner.cy(addr0, addr1);
            }

            #[pyo3(signature = (addr0, addr1, truncate = true))]
            pub fn cz(&mut self, addr0: usize, addr1: usize, truncate: bool) {
                let _ = truncate;
                self.inner.cz(addr0, addr1);
            }

            pub fn t(&mut self, addr0: usize) {
                self.inner.t(addr0);
            }

            pub fn t_adj(&mut self, addr0: usize) {
                self.inner.t_adj(addr0);
            }

            // rot1
            #[pyo3(signature = (addr0, theta, truncate = true))]
            pub fn rx(&mut self, addr0: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.rx(addr0, theta);
            }

            #[pyo3(signature = (addr0, theta, truncate = true))]
            pub fn ry(&mut self, addr0: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.ry(addr0, theta);
            }

            #[pyo3(signature = (addr0, theta, truncate = true))]
            pub fn rz(&mut self, addr0: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.rz(addr0, theta);
            }

            pub fn u3(&mut self, addr0: usize, theta: f64, phi: f64, lam: f64) {
                self.inner.u3(addr0, theta, phi, lam);
            }

            // rot2
            #[pyo3(signature = (addr0, addr1, theta, truncate = true))]
            pub fn rxx(&mut self, addr0: usize, addr1: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.rxx(addr0, addr1, theta);
            }

            #[pyo3(signature = (addr0, addr1, theta, truncate = true))]
            pub fn ryy(&mut self, addr0: usize, addr1: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.ryy(addr0, addr1, theta);
            }

            #[pyo3(signature = (addr0, addr1, theta, truncate = true))]
            pub fn rzz(&mut self, addr0: usize, addr1: usize, theta: f64, truncate: bool) {
                let _ = truncate;
                self.inner.rzz(addr0, addr1, theta);
            }

            // noise
            #[pyo3(signature = (addr0, p, truncate = true))]
            pub fn pauli_error(&mut self, addr0: usize, p: [f64; 3], truncate: bool) {
                let _ = truncate;
                self.inner.pauli_error(addr0, p);
            }

            #[pyo3(signature = (addr0, p, truncate = true))]
            pub fn depolarize(&mut self, addr0: usize, p: f64, truncate: bool) {
                let _ = truncate;
                self.inner.depolarize(addr0, p);
            }

            #[pyo3(signature = (addr0, addr1, p, truncate = true))]
            pub fn depolarize2(&mut self, addr0: usize, addr1: usize, p: f64, truncate: bool) {
                let _ = truncate;
                self.inner.depolarize2(addr0, addr1, p);
            }

            #[pyo3(signature = (addr0, addr1, p, truncate = true))]
            pub fn two_qubit_pauli_error(
                &mut self,
                addr0: usize,
                addr1: usize,
                p: [f64; 15],
                truncate: bool,
            ) {
                let _ = truncate;
                self.inner.two_qubit_pauli_error(addr0, addr1, p);
            }

            #[pyo3(signature = (addr0, p, truncate = true))]
            pub fn loss_channel(&mut self, addr0: usize, p: f64, truncate: bool) {
                let _ = truncate;
                self.inner.loss_channel(addr0, p);
            }

            #[pyo3(signature = (addr0, addr1, p, truncate = true))]
            pub fn correlated_loss_channel(
                &mut self,
                addr0: usize,
                addr1: usize,
                p: [f64; 3],
                truncate: bool,
            ) {
                let _ = truncate;
                self.inner.correlated_loss_channel(addr0, addr1, p);
            }

            #[pyo3(signature = (addr0, truncate = true))]
            pub fn reset_loss_channel(&mut self, addr0: usize, truncate: bool) {
                let _ = truncate;
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
                    ppvm_stim::sample_validated::<_, _, _, _>(
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
