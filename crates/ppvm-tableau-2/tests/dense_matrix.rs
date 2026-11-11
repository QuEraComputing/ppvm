// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Physical gates against dense matrices: row `ix_q` must equal `C† X_q C` and
//! `iz_q` must equal `C† Z_q C`, with `C` the product of textbook gate matrices.
//! No Pauli code is shared with the tableau.

use std::collections::{HashSet, VecDeque};

use num_complex::Complex64 as C;
use ppvm_tableau_2::InverseTableau;
use ppvm_traits_2::Clifford;

/// A dense `2^n × 2^n` matrix; basis index bit `q` is qubit `q`.
#[derive(Clone)]
struct Mat {
    d: usize,
    a: Vec<C>,
}

impl Mat {
    fn identity(n: usize) -> Self {
        let d = 1 << n;
        Self {
            d,
            a: (0..d * d)
                .map(|k| C::from(f64::from(u8::from(k % (d + 1) == 0))))
                .collect(),
        }
    }

    fn mul(&self, o: &Mat) -> Mat {
        let d = self.d;
        let mut a = vec![C::default(); d * d];
        for i in 0..d {
            for k in 0..d {
                let s = self.a[i * d + k];
                for j in 0..d {
                    a[i * d + j] += s * o.a[k * d + j];
                }
            }
        }
        Mat { d, a }
    }

    fn adjoint(&self) -> Mat {
        let d = self.d;
        Mat {
            d,
            a: (0..d * d)
                .map(|k| self.a[(k % d) * d + k / d].conj())
                .collect(),
        }
    }
}

type U2 = [[C; 2]; 2];

fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

/// One-qubit gate `u` on qubit `q`.
fn single(u: U2, q: usize, n: usize) -> Mat {
    let d = 1 << n;
    let mut m = Mat {
        d,
        a: vec![C::default(); d * d],
    };
    for i in 0..d {
        for j in 0..d {
            if (i ^ j) & !(1 << q) == 0 {
                m.a[i * d + j] = u[(i >> q) & 1][(j >> q) & 1];
            }
        }
    }
    m
}

/// `u` on `target`, controlled on `control` being `|1⟩`.
fn controlled(u: U2, control: usize, target: usize, n: usize) -> Mat {
    let d = 1 << n;
    let mut m = Mat {
        d,
        a: vec![C::default(); d * d],
    };
    for i in 0..d {
        for j in 0..d {
            if (i ^ j) & !(1 << target) == 0 {
                m.a[i * d + j] = if (i >> control) & 1 == 1 {
                    u[(i >> target) & 1][(j >> target) & 1]
                } else {
                    C::from(f64::from(u8::from(i == j)))
                };
            }
        }
    }
    m
}

fn gate_matrices() -> [U2; 10] {
    let (o, l, i) = (c(0.0, 0.0), c(1.0, 0.0), c(0.0, 1.0));
    let r = std::f64::consts::FRAC_1_SQRT_2;
    let h = |a: C, b: C, cc: C, dd: C| [[a * 0.5, b * 0.5], [cc * 0.5, dd * 0.5]];
    [
        [[o, l], [l, o]],                                          // X
        [[o, -i], [i, o]],                                         // Y
        [[l, o], [o, -l]],                                         // Z
        [[c(r, 0.0), c(r, 0.0)], [c(r, 0.0), c(-r, 0.0)]],         // H
        [[l, o], [o, i]],                                          // S
        [[l, o], [o, -i]],                                         // S†
        h(c(1.0, 1.0), c(1.0, -1.0), c(1.0, -1.0), c(1.0, 1.0)),   // √X
        h(c(1.0, -1.0), c(1.0, 1.0), c(1.0, 1.0), c(1.0, -1.0)),   // √X†
        h(c(1.0, 1.0), c(-1.0, -1.0), c(1.0, 1.0), c(1.0, 1.0)),   // √Y
        h(c(1.0, -1.0), c(1.0, -1.0), c(-1.0, 1.0), c(1.0, -1.0)), // √Y†
    ]
}

#[derive(Clone, Copy, Debug)]
enum Gate {
    One(usize, usize),
    Two(usize, usize, usize),
}

fn apply(t: &mut InverseTableau, gate: Gate) {
    match gate {
        Gate::One(g, q) => match g {
            0 => t.x(q),
            1 => t.y(q),
            2 => t.z(q),
            3 => t.h(q),
            4 => t.s(q),
            5 => t.s_dag(q),
            6 => t.sqrt_x(q),
            7 => t.sqrt_x_dag(q),
            8 => t.sqrt_y(q),
            _ => t.sqrt_y_dag(q),
        },
        Gate::Two(g, a, b) => match g {
            0 => t.cnot(a, b),
            1 => t.cz(a, b),
            _ => t.cy(a, b),
        },
    }
}

fn matrix(gate: Gate, n: usize) -> Mat {
    let u = gate_matrices();
    match gate {
        Gate::One(g, q) => single(u[g], q, n),
        Gate::Two(g, a, b) => controlled([u[0], u[2], u[1]][g], a, b, n),
    }
}

fn all_gates(n: usize) -> Vec<Gate> {
    let mut gates: Vec<Gate> = (0..10)
        .flat_map(|g| (0..n).map(move |q| Gate::One(g, q)))
        .collect();
    for g in 0..3 {
        for a in 0..n {
            for b in (0..n).filter(|&b| b != a) {
                gates.push(Gate::Two(g, a, b));
            }
        }
    }
    gates
}

/// All Hermitian Pauli strings on `n` qubits, indexed by `x | z << n`.
struct Paulis {
    n: usize,
    table: Vec<Mat>,
}

impl Paulis {
    fn new(n: usize) -> Self {
        let u = gate_matrices();
        let table = (0..1usize << (2 * n))
            .map(|code| {
                (0..n).fold(Mat::identity(n), |m, q| {
                    match ((code >> q) & 1, (code >> (n + q)) & 1) {
                        (0, 0) => m,
                        (1, 0) => m.mul(&single(u[0], q, n)),
                        (0, _) => m.mul(&single(u[2], q, n)),
                        _ => m.mul(&single(u[1], q, n)),
                    }
                })
            })
            .collect();
        Self { n, table }
    }

    fn x(&self, q: usize) -> &Mat {
        &self.table[1 << q]
    }

    fn z(&self, q: usize) -> &Mat {
        &self.table[1 << (self.n + q)]
    }

    /// The signed Pauli string equal to `a`, as `(x | z << n, negative)`.
    fn decompose(&self, a: &Mat) -> (usize, bool) {
        // A Pauli string maps |0⟩ to a multiple of |x⟩, so column 0 gives x.
        let x = (0..a.d)
            .find(|&r| a.a[r * a.d].norm() > 0.5)
            .expect("zero column");
        for z in 0..1usize << self.n {
            let code = x | z << self.n;
            let p = &self.table[code];
            let overlap = (0..a.d * a.d).map(|k| p.a[k].conj() * a.a[k]).sum::<C>() / a.d as f64;
            if (overlap - 1.0).norm() < 1e-9 {
                return (code, false);
            }
            if (overlap + 1.0).norm() < 1e-9 {
                return (code, true);
            }
        }
        panic!("not a signed Pauli string");
    }
}

/// Row `ix_q` (or `iz_q`) as `(x | z << n, negative)`.
fn row_code(t: &InverseTableau, q: usize, z: bool) -> (usize, bool) {
    let row = if z { t.iz(q) } else { t.ix(q) };
    let n = t.num_qubits();
    let code = (0..n)
        .map(|j| usize::from(row.x_bit(j)) << j | usize::from(row.z_bit(j)) << (n + j))
        .sum();
    (code, row.negative)
}

fn assert_rows_match(t: &InverseTableau, cliff: &Mat, paulis: &Paulis, context: &str) {
    let dag = cliff.adjoint();
    for q in 0..t.num_qubits() {
        let x = paulis.decompose(&dag.mul(paulis.x(q)).mul(cliff));
        let z = paulis.decompose(&dag.mul(paulis.z(q)).mul(cliff));
        assert_eq!(row_code(t, q, false), x, "ix_{q} after {context}");
        assert_eq!(row_code(t, q, true), z, "iz_{q} after {context}");
    }
}

/// The whole tableau packed into one integer (small `n` only).
fn key(t: &InverseTableau) -> u64 {
    let n = t.num_qubits();
    (0..n)
        .flat_map(|q| [row_code(t, q, false), row_code(t, q, true)])
        .fold(0, |k, (code, negative)| {
            (k << (2 * n + 1)) | (code as u64) << 1 | u64::from(negative)
        })
}

/// Every gate on every two-qubit Clifford. Visiting the whole group also
/// checks that the tableaux reached are exactly its 11520 elements.
#[test]
fn every_gate_on_every_two_qubit_clifford() {
    let n = 2;
    let gates = all_gates(n);
    let paulis = Paulis::new(n);
    let start = InverseTableau::new(n);
    let mut seen = HashSet::from([key(&start)]);
    let mut queue = VecDeque::from([(start, Mat::identity(n))]);
    while let Some((t, cliff)) = queue.pop_front() {
        for &gate in &gates {
            let mut next = t.clone();
            apply(&mut next, gate);
            let next_cliff = matrix(gate, n).mul(&cliff);
            assert_rows_match(&next, &next_cliff, &paulis, &format!("{gate:?}"));
            if seen.insert(key(&next)) {
                queue.push_back((next, next_cliff));
            }
        }
    }
    assert_eq!(seen.len(), 11520);
}

/// Random three-qubit circuits, covering two-qubit gates on non-adjacent qubits.
#[test]
fn random_three_qubit_circuits() {
    let n = 3;
    let gates = all_gates(n);
    let paulis = Paulis::new(n);
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    for round in 0..50 {
        let mut t = InverseTableau::new(n);
        let mut cliff = Mat::identity(n);
        for _ in 0..60 {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            let gate = gates[seed as usize % gates.len()];
            apply(&mut t, gate);
            cliff = matrix(gate, n).mul(&cliff);
        }
        assert_rows_match(&t, &cliff, &paulis, &format!("round {round}"));
    }
}
