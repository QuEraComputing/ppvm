// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// The Clifford gate set, applied in the Heisenberg picture.
/// Consumers implement the gate methods directly; aliases delegate to those methods.
pub trait Clifford {
    /// Apply Pauli `X` to one qubit.
    fn x(&mut self, qubit: usize);
    /// Apply Pauli `Y` to one qubit.
    fn y(&mut self, qubit: usize);
    /// Apply Pauli `Z` to one qubit.
    fn z(&mut self, qubit: usize);
    /// Apply Hadamard `H` to one qubit.
    fn h(&mut self, qubit: usize);
    /// Apply the phase gate `S` to one qubit.
    fn s(&mut self, qubit: usize);
    /// Apply `CNOT` to one `(control, target)` pair.
    fn cnot(&mut self, control: usize, target: usize);
    /// Apply `CZ` to one qubit pair.
    fn cz(&mut self, qubit0: usize, qubit1: usize);

    /// stim alias for [`cnot`](Clifford::cnot).
    fn cx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cnot`](Clifford::cnot).
    fn zcx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cz`](Clifford::cz).
    fn zcz(&mut self, qubit0: usize, qubit1: usize) {
        self.cz(qubit0, qubit1)
    }
    /// Apply `S†` to one qubit.
    fn s_dag(&mut self, qubit: usize);
    /// Apply `√X` to one qubit.
    fn sqrt_x(&mut self, qubit: usize);
    /// Apply `(√X)†` to one qubit.
    fn sqrt_x_dag(&mut self, qubit: usize);
    /// Apply `√Y` to one qubit.
    fn sqrt_y(&mut self, qubit: usize);
    /// Apply `(√Y)†` to one qubit.
    fn sqrt_y_dag(&mut self, qubit: usize);
    /// Apply `CY` to one `(control, target)` pair.
    fn cy(&mut self, control: usize, target: usize);
    /// stim alias for [`cy`](Clifford::cy).
    fn zcy(&mut self, control: usize, target: usize) {
        self.cy(control, target)
    }
}

/// Batched Clifford gates: apply the same gate to many qubits in one call.
pub trait CliffordBatch: Clifford {
    /// Apply Pauli `X` to every qubit in `indices`.
    fn x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.x(q);
        }
    }
    /// Apply Pauli `Y` to every qubit in `indices`.
    fn y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.y(q);
        }
    }
    /// Apply Pauli `Z` to every qubit in `indices`.
    fn z_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.z(q);
        }
    }
    /// Apply Hadamard `H` to every qubit in `indices`.
    fn h_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.h(q);
        }
    }
    /// Apply the phase gate `S` to every qubit in `indices`.
    fn s_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s(q);
        }
    }
    /// Apply `CNOT` to every `(control, target)` pair.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cnot(c, t);
        }
    }
    /// Apply `CZ` to every `(control, target)` pair.
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cz(c, t);
        }
    }
    /// apply `s†` to every qubit in `indices`.
    fn s_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s_dag(q);
        }
    }
    /// apply `√x` to every qubit in `indices`.
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x(q);
        }
    }
    /// apply `(√x)†` to every qubit in `indices`.
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x_dag(q);
        }
    }
    /// apply `√y` to every qubit in `indices`.
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y(q);
        }
    }
    /// apply `(√y)†` to every qubit in `indices`.
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y_dag(q);
        }
    }
    /// apply `cy` to every `(control, target)` pair.
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cy(c, t);
        }
    }
}
