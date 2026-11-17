// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`ApplyGateRule`]: the one trait a Pauli container implements to get [`Clifford`].
//!
//! The gate formulas live once, in `ppvm_pauli_kernels::heisenberg`. A container
//! only says how to visit its bits: it receives a whole rule and runs it in its
//! own loop, so every gate stays a single fused pass. The blanket impl below
//! maps each gate to its rule.

use ppvm_pauli_kernels::{Bits, heisenberg};

use super::Clifford;

/// A container of Pauli operators that can run a Heisenberg gate rule on them.
///
/// For each Pauli it holds, `apply1(q, rule)` reads that Pauli's bits at qubit
/// `q` as `(x, z)`, calls `rule`, writes the returned bits back, and multiplies
/// the Pauli's sign by `(−1)^flip`. With `B = bool` the container visits its
/// Paulis one at a time; with an integer `B` each bit position is a separate
/// Pauli. A container without signs ignores the flip.
///
/// With `B = bool`, `apply1` and `apply2` must call `rule` exactly once if they
/// apply the gate, and not at all if they skip it (e.g. on a lost site). Check
/// for loss before calling `rule`. Wrappers such as `Phased` take the sign flip
/// from that one call, so an extra or early call gives a wrong sign.
///
/// Implementing this trait gives [`Clifford`] through the blanket impl, so a
/// type must not also implement [`Clifford`] by hand.
pub trait ApplyGateRule {
    /// The value one rule call works on: `bool` or a bit-sliced integer.
    type B: Bits;

    /// Run a one-qubit rule on qubit `q`.
    fn apply1(&mut self, q: usize, rule: impl Fn(Self::B, Self::B) -> ([Self::B; 2], Self::B));

    /// Run a two-qubit rule on the bits `[x_a, z_a, x_b, z_b]`. For controlled
    /// gates `a` is the control.
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
