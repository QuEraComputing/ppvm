use ppvm_runtime::prelude::*;
use ppvm_runtime::strategy::CoefficientThreshold;
use pyo3::prelude::*;

#[pyclass]
pub struct TestInterface {
    inner: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>>,
}

#[pymethods]
impl TestInterface {
    #[new]
    #[pyo3(signature = (n_qubits, min_abs_coeff = 1e-10, terms = Vec::<String>::new(), coefficients = Vec::<f64>::new()))]
    pub fn new(
        n_qubits: usize,
        min_abs_coeff: f64,
        terms: Vec<String>,
        coefficients: Vec<f64>,
    ) -> Self {
        let mut ps = PauliSum::builder()
            .n_qubits(n_qubits)
            .strategy(CoefficientThreshold(min_abs_coeff))
            .capacity(n_qubits)
            .build();

        let mut coeffs = Vec::<f64>::with_capacity(terms.len());

        if coefficients.len() == 0 {
            coeffs.extend(vec![1.0; terms.len()].iter());
        } else {
            coeffs.extend(coefficients.iter());
        }

        assert_eq!(
            terms.len(),
            coeffs.len(),
            "Initial terms and coefficients need to be of same length!"
        );

        for (term, c) in terms.iter().zip(coeffs.iter()) {
            ps += (term.to_owned(), c.to_owned());
        }

        Self { inner: ps }
    }

    fn __repr__(&self) -> String {
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

    // NOTE: gate impls below
    // TODO: generate those as part of macro?

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
        self.inner.h(addr0)
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
}
