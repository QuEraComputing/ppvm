// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The blanket `Clifford` impl over `ApplyGateRule`, checked through gate identities
//! and through agreement between a `bool` and a bit-sliced `u64` container.

use ppvm_traits_2::{ApplyGateRule, Clifford};

/// One signed Pauli on `n` qubits.
#[derive(Clone, PartialEq, Debug)]
struct OnePauli {
    x: Vec<bool>,
    z: Vec<bool>,
    neg: bool,
}

impl ApplyGateRule for OnePauli {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        let ([x, z], flip) = rule(self.x[q], self.z[q]);
        (self.x[q], self.z[q]) = (x, z);
        self.neg ^= flip;
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        let ([xa, za, xb, zb], flip) = rule([self.x[a], self.z[a], self.x[b], self.z[b]]);
        (self.x[a], self.z[a], self.x[b], self.z[b]) = (xa, za, xb, zb);
        self.neg ^= flip;
    }
}

/// 64 signed Paulis on `n` qubits, bit `i` of every word belonging to Pauli `i`.
#[derive(Clone, PartialEq, Debug)]
struct Sliced {
    x: Vec<u64>,
    z: Vec<u64>,
    neg: u64,
}

impl ApplyGateRule for Sliced {
    type B = u64;

    fn apply1(&mut self, q: usize, rule: impl Fn(u64, u64) -> ([u64; 2], u64)) {
        let ([x, z], flip) = rule(self.x[q], self.z[q]);
        (self.x[q], self.z[q]) = (x, z);
        self.neg ^= flip;
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([u64; 4]) -> ([u64; 4], u64)) {
        let ([xa, za, xb, zb], flip) = rule([self.x[a], self.z[a], self.x[b], self.z[b]]);
        (self.x[a], self.z[a], self.x[b], self.z[b]) = (xa, za, xb, zb);
        self.neg ^= flip;
    }

    /// Sign-only gates skip the bit writes.
    fn apply_sign1(&mut self, q: usize, rule: impl Fn(u64, u64) -> ([u64; 2], u64)) {
        self.neg ^= rule(self.x[q], self.z[q]).1;
    }
}

impl Sliced {
    fn pauli(&self, i: usize) -> OnePauli {
        let bit = |w: u64| (w >> i) & 1 == 1;
        OnePauli {
            x: self.x.iter().map(|&w| bit(w)).collect(),
            z: self.z.iter().map(|&w| bit(w)).collect(),
            neg: bit(self.neg),
        }
    }
}

/// SplitMix64, so the tests need no dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn sliced(&mut self, n: usize) -> Sliced {
        Sliced {
            x: (0..n).map(|_| self.next()).collect(),
            z: (0..n).map(|_| self.next()).collect(),
            neg: self.next(),
        }
    }
}

const N: usize = 3;

/// Every signed Pauli on two qubits (32 of them), padded to `N` qubits.
fn all_paulis() -> Vec<OnePauli> {
    (0..32u32)
        .map(|k| {
            let bit = |i: u32| (k >> i) & 1 == 1;
            OnePauli {
                x: vec![bit(0), bit(2), false],
                z: vec![bit(1), bit(3), true],
                neg: bit(4),
            }
        })
        .collect()
}

/// `lhs` and `rhs` conjugate every Pauli the same way.
fn same_conjugation(lhs: impl Fn(&mut OnePauli), rhs: impl Fn(&mut OnePauli), what: &str) {
    for p in all_paulis() {
        let (mut a, mut b) = (p.clone(), p.clone());
        lhs(&mut a);
        rhs(&mut b);
        assert_eq!(a, b, "{what} on {p:?}");
    }
}

#[test]
fn single_qubit_identities() {
    let id = |_: &mut OnePauli| {};
    same_conjugation(|p| (p.s(0), p.s(0)).1, |p| p.z(0), "S·S = Z");
    same_conjugation(|p| (p.s(0), p.s_dag(0)).1, id, "S·S† = I");
    same_conjugation(|p| (p.s_dag(0), p.s(0)).1, id, "S†·S = I");
    same_conjugation(|p| (p.sqrt_x(0), p.sqrt_x(0)).1, |p| p.x(0), "√X·√X = X");
    same_conjugation(|p| (p.sqrt_x(0), p.sqrt_x_dag(0)).1, id, "√X·√X† = I");
    same_conjugation(|p| (p.sqrt_y(0), p.sqrt_y(0)).1, |p| p.y(0), "√Y·√Y = Y");
    same_conjugation(|p| (p.sqrt_y(0), p.sqrt_y_dag(0)).1, id, "√Y·√Y† = I");
    same_conjugation(|p| (p.h(0), p.h(0)).1, id, "H·H = I");
    same_conjugation(|p| (p.h(0), p.z(0), p.h(0)).2, |p| p.x(0), "H·Z·H = X");
    same_conjugation(|p| (p.x(0), p.y(0)).1, |p| p.z(0), "X·Y ∝ Z");
    same_conjugation(
        |p| (p.h(0), p.s(0), p.h(0)).2,
        |p| p.sqrt_x(0),
        "H·S·H ∝ √X",
    );
    // U = S·√X·S† conjugates as S first, then √X, then S†.
    same_conjugation(
        |p| (p.s(0), p.sqrt_x(0), p.s_dag(0)).2,
        |p| p.sqrt_y(0),
        "S·√X·S† ∝ √Y",
    );
}

#[test]
fn two_qubit_identities() {
    same_conjugation(
        |p| (p.h(1), p.cnot(0, 1), p.h(1)).2,
        |p| p.cz(0, 1),
        "(I⊗H)·CX·(I⊗H) = CZ",
    );
    same_conjugation(
        |p| (p.s(1), p.cnot(0, 1), p.s_dag(1)).2,
        |p| p.cy(0, 1),
        "(I⊗S)·CX·(I⊗S†) = CY",
    );
    same_conjugation(
        |p| (p.cz(0, 1), p.cz(1, 0)).1,
        |_| {},
        "CZ is symmetric and self-inverse",
    );
    same_conjugation(|p| p.cx(0, 1), |p| p.cnot(0, 1), "cx alias");
    same_conjugation(|p| p.zcy(0, 1), |p| p.cy(0, 1), "zcy alias");
}

#[test]
fn sign_only_gates_keep_bits() {
    for p in all_paulis() {
        for gate in [OnePauli::x, OnePauli::y, OnePauli::z] {
            let mut q = p.clone();
            gate(&mut q, 0);
            assert_eq!((&q.x, &q.z), (&p.x, &p.z));
        }
    }
}

/// Applies gate `g` (an index into the 13 gates) on random qubits.
fn apply_gate<C: Clifford>(c: &mut C, g: u64, a: usize, b: usize) {
    match g {
        0 => c.x(a),
        1 => c.y(a),
        2 => c.z(a),
        3 => c.h(a),
        4 => c.s(a),
        5 => c.s_dag(a),
        6 => c.sqrt_x(a),
        7 => c.sqrt_x_dag(a),
        8 => c.sqrt_y(a),
        9 => c.sqrt_y_dag(a),
        10 => c.cnot(a, b),
        11 => c.cz(a, b),
        _ => c.cy(a, b),
    }
}

#[test]
fn sliced_container_matches_single_paulis() {
    let mut rng = Rng(42);
    for _ in 0..16 {
        let mut sliced = rng.sliced(N);
        let mut singles: Vec<OnePauli> = (0..64).map(|i| sliced.pauli(i)).collect();
        for _ in 0..200 {
            let g = rng.next() % 13;
            let a = (rng.next() % N as u64) as usize;
            let b = (a + 1 + (rng.next() % (N as u64 - 1)) as usize) % N;
            apply_gate(&mut sliced, g, a, b);
            for p in &mut singles {
                apply_gate(p, g, a, b);
            }
        }
        for (i, p) in singles.iter().enumerate() {
            assert_eq!(&sliced.pauli(i), p, "Pauli {i}");
        }
    }
}
