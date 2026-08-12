// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential suite: NEW `ppvm-tableau-2` vs OLD `ppvm-tableau`.
//!
//! Both engines are driven through the single
//! [`Driver`](ppvm_conformance_2::tableau::Driver) surface, so a test can never
//! accidentally replay two different gate sequences. Everything compared here is
//! **observable algebra** — tableau rows, amplitude support and values,
//! measurement outcomes, the measurement record, expectation values — never a
//! raw hash digest (the new crate's finalization fold differs by design; the
//! hashing *contract* is tested in `tableau_hash.rs`).
//!
//! The file is organised as:
//!
//! 1. seeded-random unit differentials (construction, the full Clifford gate set
//!    row by row, the batched/fused kernels, branching gates, measurement,
//!    expectation);
//! 2. the **integration** differentials — the eight real workloads mined from the
//!    old crate's own benches and tests (MSD-85q naive and fused, the rot2
//!    brickwork golden digest, the fused-T-gate circuit, the CNOT-chain scaling
//!    sweep, the measure-all sweep, the noisy-shot average, the branch-coalesce
//!    regimes). These are the numeric acceptance bar.

use std::f64::consts::PI;

use ppvm_conformance_2::seeded_rng;
use ppvm_conformance_2::tableau::*;
use rand::RngExt;

// ===========================================================================
// helpers
// ===========================================================================

/// Every single-qubit Clifford on the `Driver` surface, as a `(name, apply)`
/// table so a test can sweep the whole gate set.
#[allow(clippy::type_complexity)]
fn single_qubit_cliffords<D: Driver>() -> Vec<(&'static str, fn(&mut D, usize))> {
    vec![
        ("x", D::x),
        ("y", D::y),
        ("z", D::z),
        ("h", D::h),
        ("s", D::s),
        ("s_dag", D::s_dag),
        ("sqrt_x", D::sqrt_x),
        ("sqrt_x_dag", D::sqrt_x_dag),
        ("sqrt_y", D::sqrt_y),
        ("sqrt_y_dag", D::sqrt_y_dag),
    ]
}

/// Every two-qubit Clifford on the `Driver` surface.
#[allow(clippy::type_complexity)]
fn two_qubit_cliffords<D: Driver>() -> Vec<(&'static str, fn(&mut D, usize, usize))> {
    vec![
        ("cnot", D::cnot),
        ("cz", D::cz),
        ("cy", D::cy),
        ("zcx", D::zcx),
        ("zcy", D::zcy),
        ("zcz", D::zcz),
    ]
}

/// A replayable gate for the randomized differentials. Written as data so the
/// identical sequence drives both engines.
#[derive(Clone, Copy, Debug)]
enum Op {
    Single(usize, usize),
    Pair(usize, usize, usize),
    T(usize),
    TDag(usize),
    Rx(usize, f64),
    Ry(usize, f64),
    Rz(usize, f64),
    Rxx(usize, usize, f64),
    Ryy(usize, usize, f64),
    Rzz(usize, usize, f64),
}

fn apply_op<D: Driver>(tab: &mut D, op: Op) {
    let singles = single_qubit_cliffords::<D>();
    let pairs = two_qubit_cliffords::<D>();
    match op {
        Op::Single(g, q) => (singles[g].1)(tab, q),
        Op::Pair(g, a, b) => (pairs[g].1)(tab, a, b),
        Op::T(q) => tab.t(q),
        Op::TDag(q) => tab.t_dag(q),
        Op::Rx(q, th) => tab.rx(q, th),
        Op::Ry(q, th) => tab.ry(q, th),
        Op::Rz(q, th) => tab.rz(q, th),
        Op::Rxx(a, b, th) => tab.rxx(a, b, th),
        Op::Ryy(a, b, th) => tab.ryy(a, b, th),
        Op::Rzz(a, b, th) => tab.rzz(a, b, th),
    }
}

/// A random circuit over the whole gate surface (Clifford + branching), emitted
/// as data.
fn random_ops(seed: u64, n_qubits: usize, len: usize, branching: bool) -> Vec<Op> {
    let mut rng = seeded_rng(seed);
    let n_single = single_qubit_cliffords::<NewWide>().len();
    let n_pair = two_qubit_cliffords::<NewWide>().len();
    let arms = if branching { 9usize } else { 2 };
    (0..len)
        .map(|_| {
            let q = rng.random_range(0..n_qubits);
            let mut b = rng.random_range(0..n_qubits);
            while b == q && n_qubits > 1 {
                b = rng.random_range(0..n_qubits);
            }
            let th: f64 = rng.random_range(-PI..PI);
            match rng.random_range(0..arms) {
                0 => Op::Single(rng.random_range(0..n_single), q),
                1 => Op::Pair(rng.random_range(0..n_pair), q, b),
                2 => Op::T(q),
                3 => Op::TDag(q),
                4 => Op::Rx(q, th),
                5 => Op::Ry(q, th),
                6 => Op::Rz(q, th),
                7 => Op::Rxx(q, b, th),
                8 => Op::Ryy(q, b, th),
                _ => Op::Rzz(q, b, th),
            }
        })
        .collect()
}

#[track_caller]
fn assert_rows_eq<A: Driver, B: Driver>(old: &A, new: &B, ctx: &str) {
    let o = old.rows();
    let n = new.rows();
    assert_eq!(o.len(), n.len(), "{ctx}: row count differs");
    for (i, (ro, rn)) in o.iter().zip(n.iter()).enumerate() {
        assert_eq!(
            ro, rn,
            "{ctx}: row {i} differs — old (x,z,phase) {ro:?} vs new {rn:?}"
        );
    }
}

#[track_caller]
fn assert_coeffs_eq<A: Driver, B: Driver>(old: &A, new: &B, tol: f64, ctx: &str) {
    let o = old.coeffs_sorted();
    let n = new.coeffs_sorted();
    assert_eq!(
        o.len(),
        n.len(),
        "{ctx}: support size differs — old {} vs new {}",
        o.len(),
        n.len()
    );
    for (a, b) in o.iter().zip(n.iter()) {
        assert_eq!(a.0, b.0, "{ctx}: support index differs");
        assert!(
            (a.1 - b.1).norm() <= tol,
            "{ctx}: coefficient at index {} differs — old {} vs new {} (tol {tol})",
            a.0,
            a.1,
            b.1
        );
    }
}

/// The full state comparison: frame rows + amplitude support/values + record +
/// loss flags.
#[track_caller]
fn assert_states_eq<A: Driver, B: Driver>(old: &A, new: &B, tol: f64, ctx: &str) {
    assert_rows_eq(old, new, ctx);
    assert_coeffs_eq(old, new, tol, ctx);
    assert_eq!(
        old.record(),
        new.record(),
        "{ctx}: measurement record differs"
    );
    assert_eq!(old.lost(), new.lost(), "{ctx}: loss flags differ");
}

// ===========================================================================
// 1. unit differentials
// ===========================================================================

#[test]
fn construction_matches() {
    for n in [1usize, 2, 5, 17, 64, 65, 85] {
        let o: OldWide = Driver::new_seeded(n, 1e-10, 3);
        let m: NewWide = Driver::new_seeded(n, 1e-10, 3);
        assert_states_eq(&o, &m, 0.0, &format!("new({n})"));
        assert_eq!(m.n_coeffs(), 1);
        assert_eq!(m.coeffs()[0].0, 0);
        assert_eq!(m.coeffs()[0].1.re, 1.0);
        assert_eq!(m.coeffs()[0].1.im, 0.0);
        assert!(m.lost().iter().all(|&l| !l));
        assert!(m.record().is_empty());
    }
}

/// Every single-qubit Clifford, applied to every qubit of a randomized frame,
/// compared row by row.
#[test]
fn single_qubit_clifford_rows_match() {
    for seed in 0..8u64 {
        let ops = random_ops(seed, 6, 24, false);
        for (gi, (name, _)) in single_qubit_cliffords::<NewWide>().iter().enumerate() {
            for q in 0..6 {
                let mut o: OldWide = Driver::new_seeded(6, 1e-10, seed);
                let mut m: NewWide = Driver::new_seeded(6, 1e-10, seed);
                for &op in &ops {
                    apply_op(&mut o, op);
                    apply_op(&mut m, op);
                }
                apply_op(&mut o, Op::Single(gi, q));
                apply_op(&mut m, Op::Single(gi, q));
                assert_rows_eq(&o, &m, &format!("seed {seed}: {name}({q})"));
            }
        }
    }
}

/// Every two-qubit Clifford on randomized frames, including the `n > 64`
/// cross-word regime.
#[test]
fn two_qubit_clifford_rows_match() {
    for &n in &[6usize, 70] {
        for seed in 0..4u64 {
            let ops = random_ops(seed, n, 20, false);
            for (gi, (name, _)) in two_qubit_cliffords::<NewWide>().iter().enumerate() {
                for &(a, b) in &[(0usize, 1usize), (1, 0), (0, n - 1), (n - 1, 0), (2, 5)] {
                    let mut o: OldWide = Driver::new_seeded(n, 1e-10, seed);
                    let mut m: NewWide = Driver::new_seeded(n, 1e-10, seed);
                    for &op in &ops {
                        apply_op(&mut o, op);
                        apply_op(&mut m, op);
                    }
                    apply_op(&mut o, Op::Pair(gi, a, b));
                    apply_op(&mut m, Op::Pair(gi, a, b));
                    assert_rows_eq(&o, &m, &format!("n={n} seed {seed}: {name}({a},{b})"));
                }
            }
        }
    }
}

/// Each batched entry point, compared old-vs-new and against its per-qubit
/// expansion, on a randomized 85-qubit frame (so word boundaries are crossed).
#[test]
fn batched_entry_points_match() {
    let n = 85;
    let idx: Vec<usize> = (0..n).step_by(3).collect();
    let pairs: Vec<(usize, usize)> = (0..17).map(|i| (i, i + 34)).collect();
    let ops = random_ops(21, n, 30, false);

    fn prep(n: usize, ops: &[Op]) -> (OldWide, NewWide) {
        let mut o: OldWide = Driver::new_seeded(n, 1e-10, 5);
        let mut m: NewWide = Driver::new_seeded(n, 1e-10, 5);
        for &op in ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        (o, m)
    }

    macro_rules! case {
        ($name:literal, $batched:expr, $expansion:expr) => {{
            let (mut o, mut m) = prep(n, &ops);
            let mut expanded = m.fork(Some(0));
            $batched(&mut o);
            $batched(&mut m);
            $expansion(&mut expanded);
            assert_rows_eq(&o, &m, concat!($name, " old-vs-new"));
            assert_rows_eq(&expanded, &m, concat!($name, " batched-vs-loop"));
        }};
    }

    case!(
        "h_many",
        |t: &mut dyn BatchDriver| t.b_h_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.h(q))
    );
    case!(
        "x_many",
        |t: &mut dyn BatchDriver| t.b_x_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.x(q))
    );
    case!(
        "s_many",
        |t: &mut dyn BatchDriver| t.b_s_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.s(q))
    );
    case!(
        "sqrt_x_many",
        |t: &mut dyn BatchDriver| t.b_sqrt_x_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.sqrt_x(q))
    );
    case!(
        "sqrt_x_dag_many",
        |t: &mut dyn BatchDriver| t.b_sqrt_x_dag_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.sqrt_x_dag(q))
    );
    case!(
        "sqrt_y_many",
        |t: &mut dyn BatchDriver| t.b_sqrt_y_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.sqrt_y(q))
    );
    case!(
        "sqrt_y_dag_many",
        |t: &mut dyn BatchDriver| t.b_sqrt_y_dag_many(&idx),
        |t: &mut NewWide| idx.iter().for_each(|&q| t.sqrt_y_dag(q))
    );
    case!(
        "cz_many",
        |t: &mut dyn BatchDriver| t.b_cz_many(&pairs),
        |t: &mut NewWide| pairs.iter().for_each(|&(a, b)| t.cz(a, b))
    );
    case!(
        "cnot_many",
        |t: &mut dyn BatchDriver| t.b_cnot_many(&pairs),
        |t: &mut NewWide| pairs.iter().for_each(|&(a, b)| t.cnot(a, b))
    );
}

/// Object-safe view of the batched surface so the macro above can call it on
/// either engine through one `dyn` reference.
trait BatchDriver {
    fn b_h_many(&mut self, qs: &[usize]);
    fn b_x_many(&mut self, qs: &[usize]);
    fn b_s_many(&mut self, qs: &[usize]);
    fn b_sqrt_x_many(&mut self, qs: &[usize]);
    fn b_sqrt_x_dag_many(&mut self, qs: &[usize]);
    fn b_sqrt_y_many(&mut self, qs: &[usize]);
    fn b_sqrt_y_dag_many(&mut self, qs: &[usize]);
    fn b_cz_many(&mut self, p: &[(usize, usize)]);
    fn b_cnot_many(&mut self, p: &[(usize, usize)]);
}
impl<T: Driver> BatchDriver for T {
    fn b_h_many(&mut self, qs: &[usize]) {
        Driver::h_many(self, qs)
    }
    fn b_x_many(&mut self, qs: &[usize]) {
        Driver::x_many(self, qs)
    }
    fn b_s_many(&mut self, qs: &[usize]) {
        Driver::s_many(self, qs)
    }
    fn b_sqrt_x_many(&mut self, qs: &[usize]) {
        Driver::sqrt_x_many(self, qs)
    }
    fn b_sqrt_x_dag_many(&mut self, qs: &[usize]) {
        Driver::sqrt_x_dag_many(self, qs)
    }
    fn b_sqrt_y_many(&mut self, qs: &[usize]) {
        Driver::sqrt_y_many(self, qs)
    }
    fn b_sqrt_y_dag_many(&mut self, qs: &[usize]) {
        Driver::sqrt_y_dag_many(self, qs)
    }
    fn b_cz_many(&mut self, p: &[(usize, usize)]) {
        Driver::cz_many(self, p)
    }
    fn b_cnot_many(&mut self, p: &[(usize, usize)]) {
        Driver::cnot_many(self, p)
    }
}

/// `cz_block` / `cz_block_pairs` — same-word, cross-word and split runs — must
/// match old *and* the per-pair `cz` expansion.
#[test]
fn cz_block_family_matches() {
    let n = 85;
    let ops = random_ops(31, n, 24, false);
    // (control_base, target_base, count): same-word, cross-word, and a run that
    // straddles a 64-bit storage-word boundary.
    for &(c, t, count) in &[
        (0usize, 17usize, 17usize),
        (0, 32, 17),
        (17, 34, 17),
        (34, 51, 17),
        (51, 68, 17),
        (60, 70, 8),
        (0, 64, 20),
    ] {
        let mut o: OldWide = Driver::new_seeded(n, 1e-10, 9);
        let mut m: NewWide = Driver::new_seeded(n, 1e-10, 9);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        let mut expanded = m.fork(Some(9));
        o.cz_block(c, t, count);
        m.cz_block(c, t, count);
        for i in 0..count {
            expanded.cz(c + i, t + i);
        }
        assert_rows_eq(&o, &m, &format!("cz_block({c},{t},{count}) old-vs-new"));
        assert_rows_eq(
            &expanded,
            &m,
            &format!("cz_block({c},{t},{count}) vs per-pair cz"),
        );
    }

    // `cz_block_pairs` (same-word entry point) on the MSD shape.
    let mut o: OldWide = Driver::new_seeded(n, 1e-10, 9);
    let mut m: NewWide = Driver::new_seeded(n, 1e-10, 9);
    for &op in &ops {
        apply_op(&mut o, op);
        apply_op(&mut m, op);
    }
    o.cz_block_pairs(0, 17, 17);
    m.cz_block_pairs(0, 17, 17);
    assert_rows_eq(&o, &m, "cz_block_pairs(0,17,17)");
}

/// Branching gates (`t`, `t_dag`, `rx`/`ry`/`rz`, `rxx`/`ryy`/`rzz`) on random
/// circuits: identical support **and** identical per-index coefficients.
#[test]
fn branching_gate_states_match() {
    for seed in 0..12u64 {
        let ops = random_ops(1000 + seed, 5, 40, true);
        let mut o: OldNarrow = Driver::new_seeded(5, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(5, 1e-10, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
            assert_rows_eq(&o, &m, &format!("seed {seed}: {op:?} rows"));
            assert_coeffs_eq(&o, &m, 1e-12, &format!("seed {seed}: {op:?} coeffs"));
        }
    }
}

/// The T-only branching path (no rot2, so no merge-order divergence): the
/// amplitude vector's **stored order** must match old element for element, not
/// merely as a set. Order is public behaviour.
#[test]
fn t_path_amplitude_order_matches() {
    for seed in 0..8u64 {
        let ops = random_ops(2000 + seed, 5, 30, true)
            .into_iter()
            .filter(|op| !matches!(op, Op::Rxx(..) | Op::Ryy(..) | Op::Rzz(..)))
            .collect::<Vec<_>>();
        let mut o: OldNarrow = Driver::new_seeded(5, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(5, 1e-10, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
            let (co, cm) = (o.coeffs(), m.coeffs());
            assert_eq!(
                co.iter().map(|e| e.0).collect::<Vec<_>>(),
                cm.iter().map(|e| e.0).collect::<Vec<_>>(),
                "seed {seed}: amplitude ORDER diverged after {op:?}"
            );
        }
    }
}

/// Seeded measurement: both sides must consume the same seed and produce the
/// identical outcome sequence *and* the identical record.
#[test]
fn measurement_outcomes_match_under_seeded_rng() {
    for seed in 0..16u64 {
        let ops = random_ops(3000 + seed, 6, 40, true);
        let mut o: OldNarrow = Driver::new_seeded(6, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(6, 1e-10, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        for q in 0..6 {
            let a = o.measure(q);
            let b = m.measure(q);
            assert_eq!(a, b, "seed {seed}: measure({q}) diverged");
            assert_states_eq(&o, &m, 1e-12, &format!("seed {seed}: after measure({q})"));
        }
        assert_eq!(o.record(), m.record(), "seed {seed}: record");
    }
}

/// Gates INTERLEAVED with measurements. The new crate keeps its measurement
/// working buffers on the tableau across `measure` calls instead of allocating a
/// fresh set per call, so the state a measurement inherits is no longer trivially
/// "whatever a constructor produced". This drives every gate class between the
/// measurements and pins the outcomes, amplitudes and record against old after
/// each step — the differential that has to hold for the reuse to be invisible.
#[test]
fn measurements_interleaved_with_gates_match_old() {
    for seed in 0..16u64 {
        let mut o: OldNarrow = Driver::new_seeded(6, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(6, 1e-10, seed);
        let ops = random_ops(7000 + seed, 6, 60, true);
        for (step, chunk) in ops.chunks(5).enumerate() {
            for &op in chunk {
                apply_op(&mut o, op);
                apply_op(&mut m, op);
            }
            let q = step % 6;
            let a = o.measure(q);
            let b = m.measure(q);
            assert_eq!(a, b, "seed {seed} step {step}: measure({q}) diverged");
            assert_states_eq(
                &o,
                &m,
                1e-12,
                &format!("seed {seed} step {step}: after measure({q})"),
            );
        }
        assert_eq!(o.record(), m.record(), "seed {seed}: record");
    }
}

/// `measure_all` and `measure_many(all)` must be observationally identical to a
/// per-qubit `measure` loop on **both** engines, and match across engines.
#[test]
fn measure_all_measure_many_and_loop_agree() {
    for seed in 0..8u64 {
        let ops = random_ops(4000 + seed, 6, 36, true);
        let mut o: OldNarrow = Driver::new_seeded(6, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(6, 1e-10, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        let all: Vec<usize> = (0..6).collect();

        let mut a1 = m.fork(Some(seed));
        let mut a2 = m.fork(Some(seed));
        let mut a3 = m.fork(Some(seed));
        let r_all = a1.measure_all();
        let r_many = a2.measure_many(&all);
        let r_loop: Vec<Option<bool>> = all.iter().map(|&q| a3.measure(q)).collect();
        assert_eq!(r_all, r_many, "seed {seed}: measure_all vs measure_many");
        assert_eq!(r_all, r_loop, "seed {seed}: measure_all vs per-qubit loop");
        assert_eq!(a1.record(), a2.record());
        assert_eq!(a1.record(), a3.record());
        assert_states_eq(&a1, &a2, 0.0, "measure_all vs measure_many state");
        assert_states_eq(&a1, &a3, 0.0, "measure_all vs loop state");

        let mut b1 = o.fork(Some(seed));
        let mut b2 = o.fork(Some(seed));
        assert_eq!(
            b1.measure_all(),
            r_all,
            "seed {seed}: old measure_all vs new"
        );
        assert_eq!(
            b2.measure_many(&all),
            r_many,
            "seed {seed}: old measure_many vs new"
        );
        assert_states_eq(&b1, &a1, 1e-12, "old vs new after measure_all");
    }
}

/// `expectation(word)` and `z_expectation(q)` on random states and random Pauli
/// words.
#[test]
fn expectation_matches() {
    for seed in 0..10u64 {
        let ops = random_ops(5000 + seed, 4, 28, true);
        let mut o: OldNarrow = Driver::new_seeded(4, 1e-10, seed);
        let mut m: NewNarrow = Driver::new_seeded(4, 1e-10, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        let mut rng = seeded_rng(seed ^ 0xabc);
        for _ in 0..24 {
            let w: String = (0..4)
                .map(|_| ['I', 'X', 'Y', 'Z'][rng.random_range(0..4usize)])
                .collect();
            let a = o.expectation_str(&w);
            let b = m.expectation_str(&w);
            assert!(
                (a - b).abs() < 1e-12,
                "seed {seed}: <{w}> old {a} vs new {b}"
            );
        }
        for q in 0..4 {
            let a = o.z_expectation(q);
            let b = m.z_expectation(q);
            assert!(
                (a - b).abs() < 1e-12,
                "seed {seed}: z_expectation({q}) old {a} vs new {b}"
            );
        }
    }
}

/// Golden expectation values from the old crate's own unit tests, reproduced on
/// the new engine (behavioural contract 15).
#[test]
fn golden_expectation_values() {
    let tol = 1e-12;
    let close =
        |a: f64, b: f64, ctx: &str| assert!((a - b).abs() < tol, "{ctx}: got {a}, expected {b}");

    let mut z0: NewNarrow = Driver::new_seeded(1, 1e-12, 0);
    close(z0.expectation_str("Z"), 1.0, "<Z> on |0>");
    close(z0.expectation_str("X"), 0.0, "<X> on |0>");
    close(z0.expectation_str("I"), 1.0, "<I> on |0>");

    let mut plus: NewNarrow = Driver::new_seeded(1, 1e-12, 0);
    plus.h(0);
    close(plus.expectation_str("X"), 1.0, "<X> on |+>");

    let mut bell: NewNarrow = Driver::new_seeded(2, 1e-12, 0);
    bell.h(0);
    bell.cnot(0, 1);
    for (w, want) in [
        ("II", 1.0),
        ("ZZ", 1.0),
        ("XX", 1.0),
        ("YY", -1.0),
        ("IZ", 0.0),
        ("ZI", 0.0),
        ("XZ", 0.0),
        ("YX", 0.0),
    ] {
        close(bell.expectation_str(w), want, &format!("Bell <{w}>"));
    }

    let mut ghz: NewNarrow = Driver::new_seeded(3, 1e-12, 0);
    ghz.h(0);
    ghz.cnot(0, 1);
    ghz.cnot(1, 2);
    for (w, want) in [
        ("III", 1.0),
        ("ZZZ", 0.0),
        ("ZIZ", 1.0),
        ("ZZI", 1.0),
        ("IZI", 0.0),
        ("XXX", 1.0),
        ("YYY", 0.0),
    ] {
        close(ghz.expectation_str(w), want, &format!("GHZ <{w}>"));
    }

    for theta in [0.0, 0.3, 1.1, PI / 2.0, 2.5] {
        let mut t: NewNarrow = Driver::new_seeded(1, 1e-12, 0);
        t.ry(0, theta);
        close(
            t.expectation_str("Z"),
            theta.cos(),
            &format!("RY({theta}) <Z>"),
        );
        close(
            t.expectation_str("X"),
            theta.sin(),
            &format!("RY({theta}) <X>"),
        );
    }

    let mut ht: NewNarrow = Driver::new_seeded(1, 1e-12, 0);
    ht.h(0);
    ht.t(0);
    let inv_sqrt2 = 1.0 / 2.0_f64.sqrt();
    close(ht.expectation_str("X"), inv_sqrt2, "H;T <X>");
    close(ht.expectation_str("Y"), inv_sqrt2, "H;T <Y>");
}

/// `Debug` must render on both engines.
#[test]
fn debug_renders() {
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 0);
    m.h(0);
    m.t(0);
    let s = format!("{m:?}");
    assert!(s.contains("GeneralizedTableau"), "Debug output: {s}");
    assert!(s.contains("coefficients"), "Debug output: {s}");
}

/// `Display` is a user-facing surface: `{}` on a tableau must render
/// **byte-for-byte** what the old crate rendered, or code printing a tableau
/// today changes output (or stops compiling). Both the bare frame and the
/// generalized tableau are compared as raw strings over a randomized circuit.
#[test]
fn display_renders_identically_to_old() {
    for seed in 0..8u64 {
        // Two-qubit rotations are excluded on purpose: old merged `rotate_2`
        // into a `std::collections::HashMap` seeded from process entropy, so
        // its amplitude ORDER — and hence the `Coefficients:` block — is
        // randomized per process and there is no old order to reproduce. That
        // is the one adjudicated divergence (see `ppvm-tableau-2`'s lib docs,
        // deferral #8); the `..._after_rotate_2` test below pins what does
        // still have to hold there.
        let mut ops = random_ops(seed, 4, 24, true);
        ops.retain(|op| !matches!(op, Op::Rxx(..) | Op::Ryy(..) | Op::Rzz(..)));
        let mut o: OldNarrow = Driver::new_seeded(4, 1e-12, seed);
        let mut m: NewNarrow = Driver::new_seeded(4, 1e-12, seed);
        for &op in &ops {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }

        // The bare Clifford frame.
        assert_eq!(
            o.tableau.to_string(),
            m.tableau.to_string(),
            "seed {seed}: Tableau Display diverged"
        );

        // The generalized tableau, verbatim: header, nested frame (including the
        // blank line after it), every `Index i: re+imi` line in order, and the
        // `Is Lost:` block.
        assert_eq!(
            o.to_string(),
            m.to_string(),
            "seed {seed}: GeneralizedTableau Display diverged"
        );
    }
}

/// After a `rotate_2` the amplitude ORDER is the one place old is
/// unreproducible (entropy-seeded `std::collections::HashMap`), so `Display`
/// can only be pinned up to that: the frame block must still match verbatim and
/// the coefficient block must match as a SET of lines.
#[test]
fn display_renders_identically_to_old_after_rotate_2() {
    for seed in 0..8u64 {
        let mut o: OldNarrow = Driver::new_seeded(3, 1e-12, seed);
        let mut m: NewNarrow = Driver::new_seeded(3, 1e-12, seed);
        for op in [Op::Single(3, 0), Op::T(0), Op::Rxx(0, 1, 0.3 * PI)] {
            apply_op(&mut o, op);
            apply_op(&mut m, op);
        }
        assert_eq!(
            o.tableau.to_string(),
            m.tableau.to_string(),
            "seed {seed}: Tableau Display diverged after rotate_2"
        );
        let (so, sn) = (o.to_string(), m.to_string());
        let mut lo: Vec<&str> = so.lines().collect();
        let mut ln: Vec<&str> = sn.lines().collect();
        lo.sort_unstable();
        ln.sort_unstable();
        assert_eq!(lo, ln, "seed {seed}: Display line set diverged");
    }
}

/// Loss flags and a non-trivial measurement history must not perturb the
/// rendering either — the `Is Lost:` block is part of old's layout.
#[test]
fn display_renders_identically_to_old_with_loss() {
    let mut o: OldNarrow = Driver::new_seeded(3, 1e-12, 11);
    let mut m: NewNarrow = Driver::new_seeded(3, 1e-12, 11);
    o.h(0);
    m.h(0);
    o.cnot(0, 1);
    m.cnot(0, 1);
    o.loss_channel(2, 1.0);
    m.loss_channel(2, 1.0);
    assert_eq!(o.tableau.to_string(), m.tableau.to_string());
    assert_eq!(o.to_string(), m.to_string());
}

// ===========================================================================
// 2. integration differentials (the numeric acceptance bar)
// ===========================================================================

/// Integration baseline #1 — MSD-85q, naive. The full 85-bit measurement
/// bitstring must be EXACTLY equal over a seed sweep, and the pre-measurement
/// amplitude support/values must match to 1e-12.
#[test]
fn integration_msd_85q_naive() {
    for seed in 0..64u64 {
        let o: OldWide = msd_state(Some(seed));
        let m: NewWide = msd_state(Some(seed));
        assert_coeffs_eq(&o, &m, 1e-12, &format!("msd seed {seed}: pre-measure"));
        assert_rows_eq(&o, &m, &format!("msd seed {seed}: pre-measure rows"));

        let so: String = msd_bitstring::<OldWide>(Some(seed));
        let sn: String = msd_bitstring::<NewWide>(Some(seed));
        assert_eq!(so.len(), 85);
        assert_eq!(so, sn, "msd seed {seed}: bitstring diverged");
    }
}

/// Integration baseline #2 — MSD-85q, fused. (a) fused == naive on the same seed
/// within the NEW crate (the old crate's own `tests/msd_batch.rs` invariant);
/// (b) fused new == fused old bit for bit; plus an exact row-by-row snapshot
/// after the Clifford portion.
#[test]
fn integration_msd_85q_fused() {
    for seed in 0..16u64 {
        let o: OldWide = msd_state_fused(Some(seed));
        let m: NewWide = msd_state_fused(Some(seed));
        assert_rows_eq(&o, &m, &format!("msd-fused seed {seed}: rows"));
        assert_coeffs_eq(&o, &m, 1e-12, &format!("msd-fused seed {seed}: coeffs"));

        // The fused Clifford portion must land on exactly the naive frame.
        let naive: NewWide = msd_state(Some(seed));
        assert_rows_eq(&naive, &m, &format!("msd seed {seed}: fused vs naive rows"));

        let fused_new = msd_bitstring_fused::<NewWide>(Some(seed));
        let naive_new = msd_bitstring::<NewWide>(Some(seed));
        assert_eq!(
            fused_new, naive_new,
            "msd seed {seed}: fused vs naive bitstring (new)"
        );
        let fused_old = msd_bitstring_fused::<OldWide>(Some(seed));
        assert_eq!(
            fused_new, fused_old,
            "msd seed {seed}: fused new vs fused old"
        );
    }
}

/// Integration baseline #3 — the rot2 brickwork. The old crate's locked
/// golden digest must be reproduced exactly, the support must be genuinely
/// branchy, and the other sizes must match old at runtime.
#[test]
fn integration_rot2_brickwork_golden_digest() {
    let n = 8;
    let tab: NewNarrow = rot2_brickwork(n, 3);
    assert!(
        tab.n_coeffs() > 8,
        "expected a branchy superposition, got {}",
        tab.n_coeffs()
    );
    let digest = measure_record_digest(&tab, n, 256);
    assert_eq!(
        digest, 0x2401_e08e_70e6_ecc8,
        "rot2 apply-path golden digest changed (got {digest:#018x})"
    );

    // ... and old reproduces the same digest on the same circuit, so the two
    // engines agree on the whole 256-seed measurement fan-out.
    let old: OldNarrow = rot2_brickwork(n, 3);
    assert_eq!(measure_record_digest(&old, n, 256), digest);
}

/// The rot2 brickwork sizes the perf gate uses, compared support-and-value with
/// old. `rotate_2`'s merge ORDER is the one adjudicated divergence (old's
/// `std::HashMap` order is process-random, so no old order exists to preserve),
/// hence the comparison is on the sorted support.
#[test]
fn integration_rot2_brickwork_support_matches() {
    for &(n, layers) in &[(8usize, 4usize), (10, 4), (12, 3)] {
        let o: OldNarrow = rot2_brickwork(n, layers);
        let m: NewNarrow = rot2_brickwork(n, layers);
        assert_rows_eq(&o, &m, &format!("rot2 n{n} l{layers}: rows"));
        assert_coeffs_eq(&o, &m, 1e-10, &format!("rot2 n{n} l{layers}: coeffs"));
        assert!(m.n_coeffs() > 8, "rot2 n{n} l{layers}: not branchy");
    }
}

/// Integration baseline #4 — the fused-T-gate circuit at 85 qubits. With
/// `fork(Some(0))` on both sides: identical support size after the T layer,
/// identical 85-entry measurement record, identical final loss flags.
#[test]
fn integration_fused_tgate_circuit() {
    for n_tgates in [8usize, 12] {
        let o_setup: OldWide = fused_tgate_setup(n_tgates);
        let m_setup: NewWide = fused_tgate_setup(n_tgates);
        assert_rows_eq(&o_setup, &m_setup, "fused-tgate setup rows");

        let mut o = o_setup.fork(Some(0));
        let mut m = m_setup.fork(Some(0));

        // Coefficients after the T layer only.
        let mut o_t = o_setup.fork(Some(0));
        let mut m_t = m_setup.fork(Some(0));
        for i in 0..n_tgates {
            o_t.t(i);
            m_t.t(i);
        }
        assert_eq!(
            o_t.n_coeffs(),
            m_t.n_coeffs(),
            "fused-{n_tgates}t: support size after the T layer"
        );
        assert_coeffs_eq(&o_t, &m_t, 1e-12, &format!("fused-{n_tgates}t: T layer"));

        let ro = fused_tgate_body(&mut o, n_tgates);
        let rm = fused_tgate_body(&mut m, n_tgates);
        assert_eq!(ro.len(), 85);
        assert_eq!(ro, rm, "fused-{n_tgates}t: measurement record");
        assert_eq!(o.lost(), m.lost(), "fused-{n_tgates}t: is_lost");
        assert_eq!(o.record(), m.record(), "fused-{n_tgates}t: record");
    }
}

/// Integration baseline #5 — the CNOT-chain scaling sweep: identical row-by-row
/// snapshot after the Clifford+T portion and identical measurement outcomes at
/// every `n`, including the `n > 64` cross-word regime.
///
/// `n ≤ 64` runs the baseline's own `usize`-index configuration; `n > 64` is
/// promoted to a `u128` index because BOTH engines compute
/// `destab_anticomm_bits |= 1 << i` over `i in 0..n`, which overflows a `usize`
/// index at `n > 64` (a pre-existing old-crate limitation the port reproduces
/// verbatim — it panics in a debug build and wraps silently in release, which is
/// why the old `tableau-scaling-128` bench only ever ran in release).
#[test]
fn integration_tableau_scaling() {
    for n in [32usize, 64] {
        let mut o: OldScaling = Driver::new_seeded(n, 1e-10, 4);
        let mut m: NewScaling = Driver::new_seeded(n, 1e-10, 4);
        scaling_prepare(&mut o);
        scaling_prepare(&mut m);
        assert_rows_eq(&o, &m, &format!("scaling n={n}: rows"));
        assert_coeffs_eq(&o, &m, 1e-12, &format!("scaling n={n}: coeffs"));

        let ro: Vec<Option<bool>> = (0..n).map(|i| o.measure(i)).collect();
        let rm: Vec<Option<bool>> = (0..n).map(|i| m.measure(i)).collect();
        assert_eq!(ro, rm, "scaling n={n}: measurement sweep");
        assert_states_eq(&o, &m, 1e-12, &format!("scaling n={n}: post-sweep"));
    }

    for n in [96usize, 128] {
        let mut o: OldWide = Driver::new_seeded(n, 1e-10, 4);
        let mut m: NewWide = Driver::new_seeded(n, 1e-10, 4);
        scaling_prepare(&mut o);
        scaling_prepare(&mut m);
        assert_rows_eq(&o, &m, &format!("scaling n={n}: rows"));
        assert_coeffs_eq(&o, &m, 1e-12, &format!("scaling n={n}: coeffs"));

        let ro: Vec<Option<bool>> = (0..n).map(|i| o.measure(i)).collect();
        let rm: Vec<Option<bool>> = (0..n).map(|i| m.measure(i)).collect();
        assert_eq!(ro, rm, "scaling n={n}: measurement sweep");
        assert_states_eq(&o, &m, 1e-12, &format!("scaling n={n}: post-sweep"));
    }
}

/// The `usize`-index shift overflow at `n > 64` is a **shared** old/new
/// behaviour, not a port divergence: both engines panic identically in a debug
/// build. Pinned here so a future "fix" on one side alone fails loudly.
#[test]
fn usize_index_overflow_is_shared_behaviour() {
    let old_panicked = std::panic::catch_unwind(|| {
        let mut o: OldScaling = Driver::new_seeded(96, 1e-10, 4);
        scaling_prepare(&mut o);
    })
    .is_err();
    let new_panicked = std::panic::catch_unwind(|| {
        let mut m: NewScaling = Driver::new_seeded(96, 1e-10, 4);
        scaling_prepare(&mut m);
    })
    .is_err();
    assert_eq!(
        old_panicked, new_panicked,
        "a `usize` index at n=96 must behave identically on both engines \
         (old panicked: {old_panicked}, new panicked: {new_panicked})"
    );
}

/// Integration baseline #6 — `measure_all` / `measure_many` on the prepared
/// 85-qubit MSD state: the three entry points agree element for element and the
/// records match, on both engines.
#[test]
fn integration_measure_all_on_msd_state() {
    let o: OldWide = msd_state(Some(3));
    let m: NewWide = msd_state(Some(3));
    let all: Vec<usize> = (0..MSD_QUBITS).collect();

    let mut m1 = m.fork(Some(11));
    let mut m2 = m.fork(Some(11));
    let mut m3 = m.fork(Some(11));
    let r_all = m1.measure_all();
    let r_many = m2.measure_many(&all);
    let r_loop: Vec<Option<bool>> = all.iter().map(|&q| m3.measure(q)).collect();
    assert_eq!(r_all, r_many);
    assert_eq!(r_all, r_loop);
    assert_eq!(m1.record(), m2.record());
    assert_eq!(m1.record(), m3.record());

    let mut o1 = o.fork(Some(11));
    let mut o2 = o.fork(Some(11));
    let mut o3 = o.fork(Some(11));
    assert_eq!(o1.measure_all(), r_all, "old measure_all vs new");
    assert_eq!(o2.measure_many(&all), r_many, "old measure_many vs new");
    let o_loop: Vec<Option<bool>> = all.iter().map(|&q| o3.measure(q)).collect();
    assert_eq!(o_loop, r_loop, "old per-qubit loop vs new");
    assert_states_eq(&o1, &m1, 1e-12, "old vs new after measure_all on MSD");
}

/// Integration baseline #7 — the 4000-shot noisy-Clifford average. (a) the
/// per-shot `⟨ZZ⟩` sequence is element-wise identical to old (which pins the RNG
/// draw order through `depolarize1` exactly); (b) the Monte-Carlo mean agrees
/// with the analytic depolarized value within `5/sqrt(N)`.
#[test]
fn integration_noisy_clifford_shots() {
    const SHOTS: u64 = 4000;
    let mut sum_new = 0.0;
    for shot in 0..SHOTS {
        let a = noisy_shot::<OldNarrow>(shot);
        let b = noisy_shot::<NewNarrow>(shot);
        assert!(
            (a - b).abs() <= 1e-12,
            "shot {shot}: <ZZ> old {a} vs new {b}"
        );
        sum_new += b;
    }
    let mean = sum_new / SHOTS as f64;
    // Three single-qubit depolarizing channels at p = 0.05 act on the ZZ
    // observable; each flips its sign with probability 2p/3.
    let f = 1.0 - 2.0 * 0.05 / 3.0;
    let expected = f * f * f;
    let tol = 5.0 / (SHOTS as f64).sqrt();
    assert!(
        (mean - expected).abs() < tol,
        "noisy-shot mean {mean} vs analytic {expected} (tol {tol})"
    );
}

/// Integration baseline #8 — the branch-coalesce regimes. For each `m = 2^j`,
/// one more `T` gate in the **doubling** regime (fresh index bit → output `2m`)
/// and in the **merge** regime (already-branched bit → output `m`): the output
/// index set and per-index coefficients must match old within 1e-9.
#[test]
fn integration_branch_coalesce_regimes() {
    for j in [2usize, 5, 8, 11] {
        // doubling: the extra T touches a qubit that has not branched yet.
        let mut o: OldWide = branch_grow(j);
        let mut m: NewWide = branch_grow(j);
        assert_eq!(o.n_coeffs(), 1 << j, "grow to 2^{j}");
        assert_eq!(m.n_coeffs(), 1 << j, "grow to 2^{j}");
        o.h(j);
        m.h(j);
        o.t(j);
        m.t(j);
        assert_eq!(o.n_coeffs(), 1 << (j + 1), "doubling regime output size");
        assert_coeffs_eq(&o, &m, 1e-9, &format!("branch doubling j={j}"));

        // merge: the extra T reuses an already-branched qubit, so every branch
        // coalesces onto an existing index.
        let mut o2: OldWide = branch_grow(j);
        let mut m2: NewWide = branch_grow(j);
        o2.t(0);
        m2.t(0);
        assert_eq!(o2.n_coeffs(), 1 << j, "merge regime output size");
        assert_coeffs_eq(&o2, &m2, 1e-9, &format!("branch merge j={j}"));
    }
}

/// Integration baseline #8, the *other* branch of the coalesce dichotomy: the
/// GENERIC `(I, u32)` fallback.
///
/// `branch_with_coefficients` takes the `u64`-packed fast path only when
/// `n_coefficients <= 0xFFFF` **and** every branch key fits in 47 bits; both
/// sequential engines otherwise fall back to a `(I, u32)` sort.
/// `integration_branch_coalesce_regimes`
/// only ever reaches `j = 11` on qubits `0..j`, so it exercises the packed path
/// exclusively and a divergence in the fallback would go unseen.
///
/// This trips the **key** predicate cheaply: branching on qubits `≥ 47` of an
/// 80-qubit `u128`-index tableau makes `idx ^ stab_anticomm_bits ≥ 2^47`, so the
/// fallback runs with a tiny support. Both the doubling and the merge regime are
/// covered, and the emitted `(index, coeff)` SEQUENCE (not just the set) is
/// compared, since the fallback has its own merge loop whose output order is
/// independently observable through the public `coefficients` field.
#[test]
fn branch_coalesce_generic_fallback_on_wide_keys() {
    // `branch_grow`'s shape, lifted onto high qubits so the keys exceed 2^47.
    fn grow_high<D: Driver>(j: usize) -> D {
        let mut tab: D = Driver::new_seeded(80, 0.0, 12345);
        for i in 0..j {
            tab.h(47 + i);
            tab.t(47 + i);
        }
        tab
    }

    for j in [1usize, 4, 8] {
        // doubling: the extra T touches a qubit that has not branched yet.
        let mut o: OldWide = grow_high(j);
        let mut m: NewWide = grow_high(j);
        assert_eq!(o.n_coeffs(), 1 << j);
        assert_eq!(m.n_coeffs(), 1 << j);
        assert!(
            o.coeffs().iter().any(|c| c.0 >= 1u128 << 47),
            "j={j}: expected a key ≥ 2^47 so the generic fallback is taken"
        );
        o.h(47 + j);
        m.h(47 + j);
        o.t(47 + j);
        m.t(47 + j);
        assert_eq!(o.n_coeffs(), 1 << (j + 1), "fallback doubling output size");
        assert_coeffs_eq(&o, &m, 1e-12, &format!("fallback doubling j={j}"));
        assert_eq!(
            o.coeffs(),
            m.coeffs(),
            "fallback doubling j={j}: emitted amplitude ORDER diverged"
        );

        // merge: the extra T reuses an already-branched qubit.
        let mut o2: OldWide = grow_high(j);
        let mut m2: NewWide = grow_high(j);
        o2.t(47);
        m2.t(47);
        assert_eq!(o2.n_coeffs(), 1 << j, "fallback merge output size");
        assert_coeffs_eq(&o2, &m2, 1e-12, &format!("fallback merge j={j}"));
        assert_eq!(
            o2.coeffs(),
            m2.coeffs(),
            "fallback merge j={j}: emitted amplitude ORDER diverged"
        );
    }
}

/// The second fallback predicate: `n_coefficients > 0xFFFF`.
///
/// The baseline's `fused-tgate` sweep tops out at 12 T gates (4096 branches) and
/// the coalesce sweep at `j = 11`, so nothing in the suite crosses the 65535
/// support cutoff. A drifted port would pass every other test and then diverge
/// on exactly the workload the cutoff exists for. `j = 16` gives a 65536-entry
/// input to the final `T`, one past the packed cutoff, while `j = 17` covers the
/// next scaling point. With Rayon enabled these also verify the large-support
/// FxHashMap path's exact public order.
///
/// Marked `#[ignore]` only because a 262144-branch doubling is slow in an
/// unoptimised test build; run it with
/// `cargo test -p ppvm-conformance-2 --release --test tableau_diff -- --ignored`.
#[test]
#[ignore = "slow in a debug build; run under --release"]
fn branch_coalesce_generic_fallback_above_the_support_cutoff() {
    for j in [16usize, 17] {
        let mut o: OldWide = branch_grow(j);
        let mut m: NewWide = branch_grow(j);
        assert_eq!(o.n_coeffs(), 1 << j);
        assert_eq!(m.n_coeffs(), 1 << j);

        o.h(j);
        m.h(j);
        o.t(j);
        m.t(j);
        assert_eq!(o.n_coeffs(), 1 << (j + 1), "doubling output size");
        assert_coeffs_eq(&o, &m, 1e-9, &format!("wide support j={j}, doubling"));
        assert_eq!(
            o.coeffs(),
            m.coeffs(),
            "wide support j={j}, doubling: emitted amplitude ORDER diverged"
        );

        let mut o2: OldWide = branch_grow(j);
        let mut m2: NewWide = branch_grow(j);
        o2.t(0);
        m2.t(0);
        assert_eq!(o2.n_coeffs(), 1 << j, "merge output size");
        assert_coeffs_eq(&o2, &m2, 1e-9, &format!("wide support j={j}, merge"));
        assert_eq!(
            o2.coeffs(),
            m2.coeffs(),
            "wide support j={j}, merge: emitted amplitude ORDER diverged"
        );
    }
}

// ===========================================================================
// 3. the remaining gate surface (not on the `Driver` shim)
// ===========================================================================

/// `u3`, `r`, the `*_many` rotation/T batches, the Pauli-error channels and the
/// basis resets — every remaining public entry point, driven with matched seeds
/// and compared state for state.
///
/// Called with fully-qualified trait paths on both sides: the two towers export
/// same-named traits, so a glob import would make every call ambiguous.
#[test]
fn remaining_gate_surface_matches() {
    use num::complex::Complex64;
    use ppvm_traits::traits as ot;
    use ppvm_traits_2 as nt;

    macro_rules! pair {
        ($seed:expr) => {{
            let o: OldNarrow = Driver::new_seeded(4, 1e-10, $seed);
            let n: NewNarrow = Driver::new_seeded(4, 1e-10, $seed);
            (o, n)
        }};
    }

    for seed in 0..8u64 {
        let th = 0.3 + seed as f64 * 0.17;

        // u3 / r
        let (mut o, mut m) = pair!(seed);
        o.h(0);
        m.h(0);
        ot::U3Gate::u3(&mut o, 0, th, th * 0.5, th * 1.5);
        nt::U3Gate::<Complex64, f64>::u3(&mut *m, 0, th, th * 0.5, th * 1.5);
        ot::RotXY::r(&mut o, 1, th, th * 0.7);
        nt::RotXY::<Complex64, f64>::r(&mut *m, 1, th, th * 0.7);
        assert_states_eq(&o, &m, 1e-12, &format!("seed {seed}: u3 + r"));

        // batched rotations and T gates
        let (mut o, mut m) = pair!(seed);
        for q in 0..4 {
            o.h(q);
            m.h(q);
        }
        let targets = [0usize, 2];
        ot::TGate::t_many(&mut o, &targets);
        nt::TGate::t_many(&mut *m, &targets);
        ot::TGate::t_dag_many(&mut o, &targets);
        nt::TGate::t_dag_many(&mut *m, &targets);
        ot::RotationOne::rx_many(&mut o, &targets, th);
        nt::RotationOne::<Complex64, f64>::rx_many(&mut *m, &targets, th);
        ot::RotationOne::ry_many(&mut o, &targets, th);
        nt::RotationOne::<Complex64, f64>::ry_many(&mut *m, &targets, th);
        ot::RotationOne::rz_many(&mut o, &targets, th);
        nt::RotationOne::<Complex64, f64>::rz_many(&mut *m, &targets, th);
        assert_states_eq(&o, &m, 1e-12, &format!("seed {seed}: batched rotations"));

        // Pauli-error channels (seeded — the RNG draw order is part of the diff)
        let (mut o, mut m) = pair!(seed);
        o.h(0);
        m.h(0);
        ot::PauliError::pauli_error(&mut o, 0, [0.1, 0.2, 0.3]);
        m.pauli_error(0, [0.1, 0.2, 0.3]);
        ot::PauliError::x_error(&mut o, 1, 0.4);
        m.x_error(1, 0.4);
        ot::PauliError::y_error(&mut o, 2, 0.4);
        m.y_error(2, 0.4);
        ot::PauliError::z_error(&mut o, 3, 0.4);
        m.z_error(3, 0.4);
        let p15 = [1.0 / 30.0; 15];
        ot::TwoQubitPauliError::two_qubit_pauli_error(&mut o, 0, 1, p15);
        m.two_qubit_pauli_error(0, 1, p15);
        ot::Depolarizing::depolarize1(&mut o, 2, 0.3);
        m.depolarize1(2, 0.3);
        ot::Depolarizing2::depolarize2(&mut o, 2, 3, 0.3);
        m.depolarize2(2, 3, 0.3);
        assert_states_eq(&o, &m, 1e-12, &format!("seed {seed}: Pauli channels"));
        assert_eq!(
            o.measure_all(),
            m.measure_all(),
            "seed {seed}: post-channel"
        );

        // basis resets
        let (mut o, mut m) = pair!(seed);
        for q in 0..4 {
            o.h(q);
            m.h(q);
            o.t(q);
            m.t(q);
        }
        ot::Reset::reset_z(&mut o, 0);
        m.reset_z(0);
        ot::Reset::reset_x(&mut o, 1);
        m.reset_x(1);
        ot::Reset::reset_y(&mut o, 2);
        m.reset_y(2);
        ot::Reset::reset_many(&mut o, &[3]);
        m.reset_many(&[3]);
        assert_states_eq(&o, &m, 1e-12, &format!("seed {seed}: basis resets"));
    }
}
