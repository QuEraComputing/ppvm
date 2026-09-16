// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::pauli::{BlanketClifford, PhaseTrack, SymplecticColumns};

/// The Clifford gate set, applied in the Heisenberg picture.
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
    /// stim alias for [`cy`](CliffordExtensions::cy).
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

impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> Clifford for T {
    #[inline]
    fn x(&mut self, q: usize) {
        self.x_phase(q);
    }

    #[inline]
    fn y(&mut self, q: usize) {
        self.y_phase(q);
    }

    #[inline]
    fn z(&mut self, q: usize) {
        self.z_phase(q);
    }

    #[inline]
    fn h(&mut self, q: usize) {
        self.flip_phase_where_xz(q);
        self.swap_xz(q);
    }

    #[inline]
    fn s(&mut self, q: usize) {
        self.s_phase(q);
        self.xor_z_from_x(q);
    }

    #[inline]
    fn cnot(&mut self, c: usize, t: usize) {
        self.cnot_phase(c, t);
        self.xor_x_col(c, t);
        self.xor_z_col(t, c);
    }

    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        self.cz_phase(a, b);
        self.cz_bits(a, b);
    }
    #[inline]
    fn s_dag(&mut self, q: usize) {
        self.s(q);
        self.z(q);
    }

    #[inline]
    fn sqrt_x(&mut self, q: usize) {
        self.h(q);
        self.s(q);
        self.h(q);
    }

    #[inline]
    fn sqrt_x_dag(&mut self, q: usize) {
        self.h(q);
        self.s_dag(q);
        self.h(q);
    }

    #[inline]
    fn sqrt_y(&mut self, q: usize) {
        self.h(q);
        self.z(q);
    }

    #[inline]
    fn sqrt_y_dag(&mut self, q: usize) {
        self.z(q);
        self.h(q);
    }

    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        self.s(target);
        self.cnot(control, target);
        self.s_dag(target);
    }
}
