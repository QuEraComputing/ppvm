// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`ApplyGateRule`]: the one trait a type implements to get [`Clifford`].
//!
//! The gate formulas live once, in `ppvm_pauli_kernels::heisenberg`. A type
//! only says how to apply a rule: it receives a whole rule and runs it in its
//! own loop, so every gate stays a single fused pass. The blanket impl below
//! maps each gate to its rule.

use ppvm_pauli_kernels::{Bits, heisenberg};

use super::Clifford;

/// A type that applies a Clifford gate given by its Heisenberg rule.
///
/// `apply1(q, rule)` applies the gate whose rule is `rule` (it maps `P` to
/// `G† P G`) at qubit `q`; `apply2` does the same for two qubits. How depends on
/// the type:
///
/// - A container of Pauli operators reads each Pauli's bits at `q` as `(x, z)`,
///   calls `rule`, writes the returned bits back, and multiplies the Pauli's sign
///   by `(−1)^flip`. With `B = bool` it visits its Paulis one at a time; with an
///   integer `B` each bit position is a separate Pauli. A container without
///   signs ignores the flip.
/// - Other types may evaluate `rule` on fixed inputs instead, e.g. an inverse
///   tableau evaluates it on the gate's generators to find which rows to
///   multiply.
///
/// A container of Pauli operators with `B = bool` must call `rule` exactly once
/// in `apply1` and `apply2` if it applies the gate, and not at all if it skips
/// it (e.g. on a lost site). Check for loss before calling `rule`. Wrappers such
/// as `Phased` take the sign flip from that one call, so an extra or early call
/// gives a wrong sign. Types that evaluate `rule` on fixed inputs are not
/// bound by this and cannot be wrapped that way.
///
/// Implementing this trait gives [`Clifford`] through the blanket impl, so a
/// type must not also implement [`Clifford`] by hand.
pub trait ApplyGateRule {
    /// The value one rule call works on: `bool` or a bit-sliced integer.
    type B: Bits;

    /// Apply a one-qubit rule at qubit `q`.
    fn apply1(&mut self, q: usize, rule: impl Fn(Self::B, Self::B) -> ([Self::B; 2], Self::B));

    /// Apply a two-qubit rule at qubits `a`, `b`; the rule's bits are
    /// `[x_a, z_a, x_b, z_b]`. For controlled gates `a` is the control.
    fn apply2(
        &mut self,
        a: usize,
        b: usize,
        rule: impl Fn([Self::B; 4]) -> ([Self::B; 4], Self::B),
    );

    /// Run a rule that leaves the bits unchanged and only flips signs (X, Y, Z).
    /// Override it to skip the bit writes.
    fn apply_sign1(
        &mut self,
        q: usize,
        rule: impl Fn(Self::B, Self::B) -> ([Self::B; 2], Self::B),
    ) {
        self.apply1(q, rule)
    }
}

impl<T: ApplyGateRule> Clifford for T {
    #[inline]
    fn x(&mut self, qubit: usize) {
        self.apply_sign1(qubit, heisenberg::x)
    }
    #[inline]
    fn y(&mut self, qubit: usize) {
        self.apply_sign1(qubit, heisenberg::y)
    }
    #[inline]
    fn z(&mut self, qubit: usize) {
        self.apply_sign1(qubit, heisenberg::z)
    }
    #[inline]
    fn h(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::h)
    }
    #[inline]
    fn s(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::s)
    }
    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        self.apply2(control, target, heisenberg::cnot)
    }
    #[inline]
    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        self.apply2(qubit0, qubit1, heisenberg::cz)
    }
    #[inline]
    fn s_dag(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::s_dag)
    }
    #[inline]
    fn sqrt_x(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::sqrt_x)
    }
    #[inline]
    fn sqrt_x_dag(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::sqrt_x_dag)
    }
    #[inline]
    fn sqrt_y(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::sqrt_y)
    }
    #[inline]
    fn sqrt_y_dag(&mut self, qubit: usize) {
        self.apply1(qubit, heisenberg::sqrt_y_dag)
    }
    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        self.apply2(control, target, heisenberg::cy)
    }
}
