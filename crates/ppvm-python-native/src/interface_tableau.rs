use bnum::types::{U256, U512, U1024, U2048, U4096};
use paste::paste;
use ppvm_runtime::prelude::*;
use ppvm_runtime::tableau::GeneralizedTableau;
use pyo3::prelude::*;

macro_rules! create_interface {
    ($name: ident, $type: ident, $indexType: ident) => {
        #[pyclass]
        pub struct $name {
            inner: GeneralizedTableau<$type, $indexType>,
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10))]
            pub fn new(n_qubits: usize, min_abs_coeff: f64) -> Self {
                let tab: GeneralizedTableau<$type, $indexType> =
                    GeneralizedTableau::new(n_qubits, min_abs_coeff);

                Self { inner: tab }
            }

            fn __repr__(&self) -> String {
                // TODO: expose some more details e.g. for debugging
                format!("{}", self.inner)
            }

            fn __str__(&self) -> String {
                self.inner.to_string()
            }

            pub fn measure(&mut self, addr0: usize) -> bool {
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

            pub fn cnot(&mut self, addr0: usize, addr1: usize) {
                self.inner.cnot(addr0, addr1);
            }

            pub fn cz(&mut self, addr0: usize, addr1: usize) {
                self.inner.cz(addr0, addr1);
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

            pub fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [f64; 15]) {
                self.inner.two_qubit_pauli_error(addr0, addr1, p);
            }

            // some python niceties

            fn __copy__(&self) -> Self {
                Self {
                    inner: self.inner.clone(),
                }
            }

            // NOTE: PartialEq not implemented for Tableau
            // fn __richcmp__(
            //     &self,
            //     other: PyRef<$name>,
            //     op: pyo3::basic::CompareOp,
            // ) -> PyResult<bool> {
            //     match op {
            //         pyo3::basic::CompareOp::Eq => Ok(self.inner == other.inner),
            //         pyo3::basic::CompareOp::Ne => Ok(self.inner != other.inner),
            //         _ => Err(pyo3::exceptions::PyNotImplementedError::new_err(
            //             "Only equality and inequality comparisons are supported for PauliSum.",
            //         )),
            //     }
            // }

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
            type [<$name$n>] = config::indexmap::ByteFxHashF64<{(2 as usize).pow($n)}>;
            create_interface!([<GeneralizedTableau$name$n>], [<$name$n>], $indexType);
        )*
    }
    };
}

// up to 64 qubits
create_interface_range!(IndexMapFxHash, usize, 0, 1, 2, 3, 4, 5, 6, 7, 8);

// 64 - 128 qubits
create_interface_range!(IndexMapFxHash, u128, 9, 10, 11, 12, 13, 14, 15, 16);

// 128 - 256 qubits
create_interface_range!(
    IndexMapFxHash,
    U256,
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
