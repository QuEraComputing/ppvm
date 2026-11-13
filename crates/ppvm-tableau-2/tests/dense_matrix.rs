// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Gates against dense matrices: row `ix_q` must equal `C† X_q C` and
//! `iz_q` must equal `C† Z_q C`, with `C` the product of textbook gate matrices.
//! Measurements are checked against a state vector. No Pauli code is shared
//! with the tableau.

use std::collections::{HashSet, VecDeque};

use num_complex::Complex64 as C;
use ppvm_tableau_2::{ActiveDormant, DenseState, GeneralizedTableau, InverseTableau};
use ppvm_traits_2::{Clifford, Pauli, RotationOne, RotationTwo, TGate};

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

fn apply(t: &mut impl Clifford, gate: Gate) {
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

type Rule1 = fn(bool, bool) -> ([bool; 2], bool);
type Rule2 = fn([bool; 4]) -> ([bool; 4], bool);

/// `gate` as a change of basis on virtual qubits, `C ← C·V`.
fn change_basis(t: &mut InverseTableau, gate: Gate) {
    use ppvm_pauli_kernels::heisenberg as rules;
    let one: [Rule1; 10] = [
        rules::x,
        rules::y,
        rules::z,
        rules::h,
        rules::s,
        rules::s_dag,
        rules::sqrt_x,
        rules::sqrt_x_dag,
        rules::sqrt_y,
        rules::sqrt_y_dag,
    ];
    let two: [Rule2; 3] = [rules::cnot, rules::cz, rules::cy];
    match gate {
        Gate::One(g, v) => t.change_basis1(v, one[g]),
        Gate::Two(g, a, b) => t.change_basis2(a, b, two[g]),
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

/// Every gate, physical and as a change of basis, on every two-qubit
/// Clifford. Visiting the whole group also
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
            let mut rebased = t.clone();
            change_basis(&mut rebased, gate);
            let rebased_cliff = cliff.mul(&matrix(gate, n));
            assert_rows_match(
                &rebased,
                &rebased_cliff,
                &paulis,
                &format!("basis {gate:?}"),
            );
            if seen.insert(key(&next)) {
                queue.push_back((next, next_cliff));
            }
        }
    }
    assert_eq!(seen.len(), 11520);
}

/// Random three-qubit circuits mixing physical gates and changes of basis,
/// covering two-qubit gates on non-adjacent qubits.
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
            if seed >> 63 == 0 {
                apply(&mut t, gate);
                cliff = matrix(gate, n).mul(&cliff);
            } else {
                change_basis(&mut t, gate);
                cliff = cliff.mul(&matrix(gate, n));
            }
        }
        assert_rows_match(&t, &cliff, &paulis, &format!("round {round}"));
    }
}

/// Row `ix_q` (or `iz_q`) as a signed Pauli matrix.
fn row_matrix(t: &InverseTableau, q: usize, z: bool, paulis: &Paulis) -> Mat {
    let (code, negative) = row_code(t, q, z);
    let mut m = paulis.table[code].clone();
    if negative {
        m.a.iter_mut().for_each(|e| *e = -*e);
    }
    m
}

/// `⟨0…0| C⁻¹ P C |0…0⟩` for the Pauli string `code`, from the rows alone.
fn tableau_expectation(t: &InverseTableau, code: usize, paulis: &Paulis) -> C {
    conjugated(t, code, paulis).a[0]
}

/// `C⁻¹ P C` for the Pauli string `code`, from the rows alone:
/// `P = i^#Y · ∏ X_q^x Z_q^z`, so `C⁻¹ P C = i^#Y · ∏ ix_q^x iz_q^z`.
fn conjugated(t: &InverseTableau, code: usize, paulis: &Paulis) -> Mat {
    let n = t.num_qubits();
    let mut m = Mat::identity(n);
    let mut phase = c(1.0, 0.0);
    for q in 0..n {
        let (x, z) = ((code >> q) & 1 == 1, (code >> (n + q)) & 1 == 1);
        if x {
            m = m.mul(&row_matrix(t, q, false, paulis));
        }
        if z {
            m = m.mul(&row_matrix(t, q, true, paulis));
        }
        if x && z {
            phase *= c(0.0, 1.0);
        }
    }
    m.a.iter_mut().for_each(|e| *e *= phase);
    m
}

fn apply_to_state(m: &Mat, v: &[C]) -> Vec<C> {
    (0..m.d)
        .map(|i| (0..m.d).map(|j| m.a[i * m.d + j] * v[j]).sum())
        .collect()
}

fn state_expectation(v: &[C], p: &Mat) -> C {
    let pv = apply_to_state(p, v);
    v.iter().zip(&pv).map(|(a, b)| a.conj() * b).sum()
}

/// Probability that qubit `q` is `|1⟩`.
fn probability_one(v: &[C], q: usize) -> f64 {
    v.iter()
        .enumerate()
        .filter(|(i, _)| (i >> q) & 1 == 1)
        .map(|(_, a)| a.norm_sqr())
        .sum()
}

/// Project onto outcome `m` of qubit `q` and renormalize.
fn project(v: &mut [C], q: usize, m: bool) {
    v.iter_mut()
        .enumerate()
        .filter(|(i, _)| ((i >> q) & 1 == 1) != m)
        .for_each(|(_, a)| *a = C::default());
    let norm = v.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
    v.iter_mut().for_each(|a| *a /= norm);
}

/// Measure `q` on both the tableau and the state vector, forcing `forced` when
/// the outcome is random. The tableau must call the coin exactly when the
/// state vector says the outcome is random.
fn measure_both(t: &mut InverseTableau, v: &mut [C], q: usize, forced: bool) -> bool {
    let p1 = probability_one(v, q);
    let random = p1 > 1e-9 && p1 < 1.0 - 1e-9;
    if random {
        assert!(
            (p1 - 0.5).abs() < 1e-9,
            "random outcome must be 50/50, got {p1}"
        );
    }
    let mut called = false;
    let m = t.measure(q, || {
        called = true;
        forced
    });
    assert_eq!(called, random, "coin called iff the outcome is random");
    if !random {
        assert_eq!(m, p1 > 0.5, "deterministic outcome");
    }
    project(v, q, m);
    m
}

/// Random three-qubit circuits with measurements and resets: after every step,
/// all 64 Pauli expectations from the tableau match the state vector.
#[test]
fn measurements_match_the_state_vector() {
    let n = 3;
    let gates = all_gates(n);
    let paulis = Paulis::new(n);
    let x = gate_matrices()[0];
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    for round in 0..30 {
        let mut t = InverseTableau::new(n);
        let mut v = vec![C::default(); 1 << n];
        v[0] = c(1.0, 0.0);
        for step in 0..60 {
            let (q, forced) = (next() % n, next() % 2 == 1);
            match next() % 10 {
                0..7 => {
                    let gate = gates[next() % gates.len()];
                    apply(&mut t, gate);
                    v = apply_to_state(&matrix(gate, n), &v);
                }
                7..9 => {
                    measure_both(&mut t, &mut v, q, forced);
                }
                _ => {
                    let m = measure_both(&mut t, &mut v, q, forced);
                    let mut t_reset = t.clone();
                    t_reset.reset(q, || unreachable!("qubit was just measured"));
                    if m {
                        t.x(q);
                        v = apply_to_state(&single(x, q, n), &v);
                    }
                    assert_eq!(t, t_reset, "reset = measure + conditional X");
                }
            }
            for code in 0..1usize << (2 * n) {
                let (from_rows, from_state) = (
                    tableau_expectation(&t, code, &paulis),
                    state_expectation(&v, &paulis.table[code]),
                );
                assert!(
                    (from_rows - from_state).norm() < 1e-9,
                    "round {round} step {step}: Pauli {code:#b}: tableau {from_rows}, state {from_state}"
                );
            }
        }
    }
}

type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

/// `|χ⟩ = |φ⟩ ⊗ |0…0⟩` as a vector over all virtual qubits.
fn virtual_state(s: &Sim) -> Vec<C> {
    let (active, phi) = (s.state().active(), s.state().phi().amplitudes());
    let mut chi = vec![C::default(); 1 << s.num_qubits()];
    for (alpha, &amp) in phi.iter().enumerate() {
        let i: usize = (0..active.len())
            .map(|j| ((alpha >> j) & 1) << active[j])
            .sum();
        chi[i] = amp;
    }
    chi
}

/// Code (`x | z << n`) of the Pauli `axis` on qubit `q`.
fn pauli_code(axis: Pauli, q: usize, n: usize) -> usize {
    let (x, z) = match axis {
        Pauli::I => (0, 0),
        Pauli::X => (1, 0),
        Pauli::Y => (1, 1),
        Pauli::Z => (0, 1),
    };
    x << q | z << (n + q)
}

/// `exp(−i θ/2 P) = cos(θ/2)·I − i sin(θ/2)·P`.
fn rotation(p: &Mat, theta: f64) -> Mat {
    let (s, co) = (theta / 2.0).sin_cos();
    let mut m = p.clone();
    for (k, e) in m.a.iter_mut().enumerate() {
        let identity = f64::from(u8::from(k % (p.d + 1) == 0));
        *e = c(co * identity, 0.0) + c(0.0, -s) * *e;
    }
    m
}

/// Random three-qubit circuits of Cliffords, `T`, `T†` and one- and two-qubit
/// Pauli rotations: after every gate, all 64 Pauli expectations of
/// `C |χ⟩` match the state vector.
#[test]
fn rotations_match_the_state_vector() {
    let n = 3;
    let gates = all_gates(n);
    let paulis = Paulis::new(n);
    let axes = [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];
    let w = C::from_polar(1.0, std::f64::consts::FRAC_PI_4);
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize
    };
    let mut max_axes = 0;
    for round in 0..40 {
        let mut s = Sim::new_with_seed(n, round);
        let mut v = vec![C::default(); 1 << n];
        v[0] = c(1.0, 0.0);
        for step in 0..40 {
            let (a, theta) = (next() % n, (next() % 1000) as f64 * 0.00628);
            let b = (a + 1 + next() % (n - 1)) % n;
            let (pa, pb) = (axes[next() % 4], axes[next() % 4]);
            let m = match next() % 8 {
                0..3 => {
                    let gate = gates[next() % gates.len()];
                    apply(&mut s, gate);
                    matrix(gate, n)
                }
                3 => {
                    s.t(a);
                    single([[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), w]], a, n)
                }
                4 => {
                    s.t_dag(a);
                    single([[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), w.conj()]], a, n)
                }
                5 | 6 => {
                    s.rotate_1(pa, a, theta);
                    rotation(&paulis.table[pauli_code(pa, a, n)], theta)
                }
                _ => {
                    s.rotate_2(pa, pb, a, b, theta);
                    let code = pauli_code(pa, a, n) | pauli_code(pb, b, n);
                    rotation(&paulis.table[code], theta)
                }
            };
            v = apply_to_state(&m, &v);
            max_axes = max_axes.max(s.state().active().len());
            let chi = virtual_state(&s);
            for code in 0..1usize << (2 * n) {
                let (from_sim, from_state) = (
                    state_expectation(&chi, &conjugated(s.tableau(), code, &paulis)),
                    state_expectation(&v, &paulis.table[code]),
                );
                assert!(
                    (from_sim - from_state).norm() < 1e-9,
                    "round {round} step {step}: Pauli {code:#b}: simulator {from_sim}, state {from_state}"
                );
            }
        }
    }
    assert_eq!(max_axes, n, "the circuits activate every qubit");
}
