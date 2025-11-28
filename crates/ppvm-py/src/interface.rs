use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
use pyo3::prelude::*;

// adapted from https://pyo3.rs/v0.27.1/class.html#no-generic-parameters
macro_rules! create_interface {
    ($name: ident, $type: ident) => {
        #[pyclass]
        pub struct $name {
            inner: PauliSum<$type>,
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, max_pauli_weight = usize::MAX, terms = Vec::<String>::new(), coefficients = Vec::<f64>::new()))]
            pub fn new(
                n_qubits: usize,
                min_abs_coeff: f64,
                max_pauli_weight: usize,
                terms: Vec<String>,
                coefficients: Vec<f64>
            ) -> Self {

                // TODO: this is not ideal since we could skip one of the strategies completely; need to look into
                // how we can do this in the macro here
                let strat = CombinedStrategy(CoefficientThreshold(min_abs_coeff), MaxPauliWeight(max_pauli_weight));
                let mut ps = PauliSum::builder()
                    .n_qubits(n_qubits)
                    .strategy(strat)
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
        }

    };
}

// TODO: use macro to loop over these? Let's check performance first
type IndexMapFxHash10 =
    config::indexmap::ByteFxHashF64<10, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
type IndexMapFxHash100 =
    config::indexmap::ByteFxHashF64<100, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
type IndexMapFxHash500 =
    config::indexmap::ByteFxHashF64<500, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;
type IndexMapFxHash1000 =
    config::indexmap::ByteFxHashF64<1000, CombinedStrategy<CoefficientThreshold, MaxPauliWeight>>;

// NOTE: let's hope no one wants to simulate more than 8k qubits
create_interface!(PauliSumIndexMapFxHash10, IndexMapFxHash10);
create_interface!(PauliSumIndexMapFxHash100, IndexMapFxHash100);
create_interface!(PauliSumIndexMapFxHash500, IndexMapFxHash500);
create_interface!(PauliSumIndexMapFxHash1000, IndexMapFxHash1000);
