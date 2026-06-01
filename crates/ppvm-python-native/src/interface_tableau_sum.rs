use bnum::types::{U256, U512, U1024, U2048};
use paste::paste;
use ppvm_runtime::prelude::*;
use ppvm_tableau_sum::prelude::*;
use ppvm_tableau_sum::sampler::Sampler;
use pyo3::prelude::*;

use crate::interface_tableau::measurement_to_u8;

macro_rules! create_sum_interface {
    ($tab_name: ident, $sampler_name: ident, $type: ident, $indexType: ident) => {
        #[pyclass]
        pub struct $tab_name {
            inner: GeneralizedTableauSum<$type, $indexType>,
        }

        #[pymethods]
        impl $tab_name {
            #[new]
            #[pyo3(signature = (
                n_qubits,
                min_abs_coeff = 1e-10,
                sum_cutoff = 1e-8,
                seed = None,
            ))]
            pub fn new(
                n_qubits: usize,
                min_abs_coeff: f64,
                sum_cutoff: f64,
                seed: Option<u64>,
            ) -> Self {
                let inner: GeneralizedTableauSum<$type, $indexType> = match seed {
                    Some(s) => GeneralizedTableauSum::new_with_seed(
                        n_qubits,
                        min_abs_coeff,
                        sum_cutoff,
                        s,
                    ),
                    None => GeneralizedTableauSum::new(n_qubits, min_abs_coeff, sum_cutoff),
                };
                Self { inner }
            }

            /// Number of branches currently in the sum.
            pub fn __len__(&self) -> usize {
                self.inner.len()
            }

            pub fn len(&self) -> usize {
                self.inner.len()
            }

            /// Mid-circuit measurement probabilities `(p_zero, p_one, p_lost)`.
            ///
            /// Unlike `GeneralizedTableau.measure` (which samples and returns a
            /// single outcome), the sum-form measurement branches each entry
            /// into its three Z-basis outcomes and returns the aggregated
            /// probabilities — use this for analytic measurement statistics,
            /// and a `Sampler` from `.sampler()` for stochastic samples.
            pub fn measure(&mut self, addr0: usize) -> (f64, f64, f64) {
                self.inner.measure(addr0)
            }

            // Clifford
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

            // Single-qubit rotations
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

            // Two-qubit rotations
            pub fn rxx(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rxx(addr0, addr1, theta);
            }

            pub fn ryy(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.ryy(addr0, addr1, theta);
            }

            pub fn rzz(&mut self, addr0: usize, addr1: usize, theta: f64) {
                self.inner.rzz(addr0, addr1, theta);
            }

            // Noise
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

            /// Force a truncation pass.
            ///
            /// Truncation is already applied automatically inside every noise
            /// op; this is exposed for users who want to prune after a long
            /// run of unitary gates.
            pub fn truncate(&mut self) {
                self.inner.truncate();
            }

            /// Compile a `Sampler` snapshotting the current state.
            ///
            /// The returned sampler holds its own RNG and a sorted copy of the
            /// sum's branches; further gates / noise on this tableau do not
            /// affect it. Call `.sample_shots(n)` on the result to draw shots.
            /// Two samplers compiled in a row use independent RNG sequences.
            pub fn sampler(&mut self) -> $sampler_name {
                $sampler_name {
                    inner: self.inner.sampler(),
                }
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

        #[pyclass]
        pub struct $sampler_name {
            inner: Sampler<$type, $indexType>,
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

macro_rules! create_sum_interface_range {
    ($indexType: ident, $( $n: expr),* ) => {
        paste! {
        $(
            type [<SumConfig$n>] = config::fx64hash::Byte8F64<$n>;
            create_sum_interface!(
                [<GeneralizedTableauSum$n>],
                [<TableauSumSampler$n>],
                [<SumConfig$n>],
                $indexType
            );
        )*
        }
    };
}

// up to 64 qubits
create_sum_interface_range!(usize, 1);

// 64 - 128 qubits
create_sum_interface_range!(u128, 2);

// 128 - 256 qubits
create_sum_interface_range!(U256, 3, 4);

create_sum_interface_range!(U512, 5, 6, 7, 8);

create_sum_interface_range!(U1024, 9, 10, 11, 12, 13, 14, 15, 16);

create_sum_interface_range!(
    U2048, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32
);
