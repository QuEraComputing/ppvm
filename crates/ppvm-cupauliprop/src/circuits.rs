// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Trotter circuits mirroring `ppvm-benchmarks/trotter-benchmarks`
//! (`python-benchmarks/ppvm/test_trotter.py` and `test_heisenberg.py`).
//!
//! Both walk the forward circuit in reverse (Heisenberg picture): ppvm gates
//! already conjugate by the adjoint, so only the operation order is flipped.

/// The gate set the benchmark circuits need, implemented by every backend.
///
/// Every gate is followed by truncation, like `ppvm.PauliSum` with the
/// default `truncate=True`.
pub trait Propagator {
    fn rx(&mut self, q: usize, theta: f64);
    fn rxx(&mut self, a: usize, b: usize, theta: f64);
    fn ryy(&mut self, a: usize, b: usize, theta: f64);
    fn rzz(&mut self, a: usize, b: usize, theta: f64);
    /// Single-qubit Pauli channel with probabilities `[p_x, p_y, p_z]`.
    fn pauli_error(&mut self, q: usize, p: [f64; 3]);
    /// `Tr(O |0…0⟩⟨0…0|)`.
    fn overlap_with_zero(&self) -> f64;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Dense Pauli strings of the benchmark observable `Σ_i Z_i`.
pub fn sum_z_terms(n_qubits: usize) -> Vec<String> {
    (0..n_qubits)
        .map(|i| {
            (0..n_qubits)
                .map(|j| if i == j { 'Z' } else { 'I' })
                .collect()
        })
        .collect()
}

/// 1D transverse-field Ising Trotter circuit (`test_trotter.py::trotter`).
///
/// Forward step: `[RX, noise]` on every site, then `[RZZ, noise, noise]` on
/// every bond `(i, i+1)`.
pub fn tfim<P: Propagator>(
    state: &mut P,
    n_qubits: usize,
    n_steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
) {
    for _ in 0..n_steps {
        for i in (0..n_qubits - 1).rev() {
            state.pauli_error(i + 1, noise);
            state.pauli_error(i, noise);
            state.rzz(i, i + 1, theta_zz);
        }
        for i in (0..n_qubits).rev() {
            state.pauli_error(i, noise);
            state.rx(i, theta_x);
        }
    }
}

/// Bonds of an `l × l` open-boundary square lattice, in the row-major order
/// of `test_heisenberg.py::lattice_edges`.
pub fn lattice_edges(l: usize) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for i in 0..l {
        for j in 0..l {
            let s = i * l + j;
            if j < l - 1 {
                edges.push((s, s + 1));
            }
            if i < l - 1 {
                edges.push((s, s + l));
            }
        }
    }
    edges
}

/// Heisenberg-model Trotter circuit (`test_heisenberg.py::heisenberg`).
///
/// Forward bond step: `[Rxx, N, N, Ryy, N, N, Rzz, N, N]`.
pub fn heisenberg<P: Propagator>(
    state: &mut P,
    edges: &[(usize, usize)],
    n_steps: usize,
    theta: f64,
    noise: [f64; 3],
) {
    for _ in 0..n_steps {
        for &(a, b) in edges.iter().rev() {
            state.pauli_error(b, noise);
            state.pauli_error(a, noise);
            state.rzz(a, b, theta);
            state.pauli_error(b, noise);
            state.pauli_error(a, noise);
            state.ryy(a, b, theta);
            state.pauli_error(b, noise);
            state.pauli_error(a, noise);
            state.rxx(a, b, theta);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Records the operation sequence instead of simulating it.
    #[derive(Default)]
    struct Recorder(Vec<String>);

    impl Propagator for Recorder {
        fn rx(&mut self, q: usize, _: f64) {
            self.0.push(format!("rx{q}"));
        }
        fn rxx(&mut self, a: usize, b: usize, _: f64) {
            self.0.push(format!("rxx{a},{b}"));
        }
        fn ryy(&mut self, a: usize, b: usize, _: f64) {
            self.0.push(format!("ryy{a},{b}"));
        }
        fn rzz(&mut self, a: usize, b: usize, _: f64) {
            self.0.push(format!("rzz{a},{b}"));
        }
        fn pauli_error(&mut self, q: usize, _: [f64; 3]) {
            self.0.push(format!("pe{q}"));
        }
        fn overlap_with_zero(&self) -> f64 {
            0.0
        }
        fn len(&self) -> usize {
            self.0.len()
        }
    }

    #[test]
    fn tfim_step_order() {
        let mut r = Recorder::default();
        tfim(&mut r, 3, 1, 0.1, 0.1, [0.0; 3]);
        let expected = [
            "pe2", "pe1", "rzz1,2", "pe1", "pe0", "rzz0,1", "pe2", "rx2", "pe1", "rx1", "pe0",
            "rx0",
        ];
        assert_eq!(r.0, expected);
    }

    #[test]
    fn tfim_op_count() {
        let mut r = Recorder::default();
        tfim(&mut r, 128, 20, 0.1, 0.1, [0.0; 3]);
        assert_eq!(r.0.len(), 20 * (3 * 127 + 2 * 128));
    }

    #[test]
    fn lattice_edges_match_python() {
        assert_eq!(lattice_edges(2), [(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(lattice_edges(6).len(), 60);
    }

    #[test]
    fn heisenberg_bond_order() {
        let mut r = Recorder::default();
        heisenberg(&mut r, &[(0, 1)], 1, 0.1, [0.0; 3]);
        let expected = [
            "pe1", "pe0", "rzz0,1", "pe1", "pe0", "ryy0,1", "pe1", "pe0", "rxx0,1",
        ];
        assert_eq!(r.0, expected);
    }

    #[test]
    fn sum_z_terms_layout() {
        assert_eq!(sum_z_terms(3), ["ZII", "IZI", "IIZ"]);
    }
}
