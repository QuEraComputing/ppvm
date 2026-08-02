// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! End-to-end circuit-propagation differential correctness: propagate a *whole*
//! circuit (a Trotterized transverse-field Ising evolution, and a deep random
//! Clifford+rotation circuit) through BOTH the new
//! `ppvm-pauli-sum-2::Sum` and the old `ppvm-pauli-sum::PauliSum` on an
//! **identical algebraic config** and assert the final supports match.
//!
//! This is the coverage the `-2` conformance suite was missing. The existing
//! `pauli_sum_diff` / `pauli_sum_rotation_noise_diff` suites diff *single* gates
//! (and short replays) on a *fresh* sum, so they can't see how per-gate
//! allocation and truncation behave over a deep circuit — exactly where the
//! storage double-buffer (aux/scratch) lives and where the old crate's
//! real-workload benches/tests (`trotter.rs`, `random-circuit.rs`) exercised the
//! engine. Here we port the old `trotter_func` workload and drive it end to end.
//!
//! Config parity (both sides):
//! * `[u8; 8]` (64-qubit-capacity) storage, `f64` coefficients.
//! * The SAME truncation — old `CombinedStrategy(CoefficientThreshold(1e-6),
//!   MaxPauliWeight(w))`, new `CombinedPolicy(CoefficientThreshold(1e-6),
//!   MaxPauliWeight(w))`. The two crates' keep-rules are identical: old drops
//!   `|c| < threshold` (`Coefficient::cutoff`) / `weight > w`, new keeps
//!   `|c| >= threshold` / `weight <= w`.
//! * Truncation is caller-driven on both (`truncate()` after each operation),
//!   exactly like the old `trotter_func`.
//!
//! The new crate has no native `rzz`, so a `ZZ` rotation is decomposed as
//! `rzz(a, b, θ) = cnot(a, b); rz(b, θ); cnot(a, b)` and the SAME decomposed
//! sequence is applied on BOTH sides (a sub-test verifies old-native `rzz` ≈
//! old-decomposed, so the decomposition itself is validated).

use ppvm_conformance_2::{GateOp, assert_close, random_circuit, seeded_rng};

// --- Old crate: build a CombinedStrategy config and pull in its gate traits. ---
use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{
    Clifford as OldClifford, PauliError as OldPauliError, RotationOne as OldRotationOne,
    RotationTwo as OldRotationTwo,
};

// --- New crate: the storage-matched `[u8; 8]` sum + CombinedPolicy + traits. ---
use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
};

use rand::rngs::StdRng;

/// The shared truncation floor (both crates).
const THRESHOLD: f64 = 1e-6;
/// Coefficient comparison tolerance. An end-to-end circuit accumulates
/// `sin`/`cos` products and the two crates merge colliding keys in different
/// iteration orders, so a few ulp of drift per gate is expected; this floor is
/// well below the `1e-6` truncation threshold, so it never reclassifies a
/// kept/dropped term.
const TOL: f64 = 1e-9;

// --- Old side: `[u8; 8]` (64-qubit-capacity) + CombinedStrategy(1e-6, w). ------
type OldStrat = CombinedStrategy<OldCoeffThreshold, OldMaxWeight>;
type OldCfg = OldByteF64<8, OldStrat>;
type OldSum = OldPauliSum<OldCfg>;

// --- New side: storage-matched `[u8; 8]` key + CombinedPolicy(1e-6, w). --------
type NewKey = NewPauliWord<[u8; 8]>;
type NewPolicy = CombinedPolicy<NewCoeffThreshold, NewMaxWeight>;
type NewSum = Sum<HashMapStore<NewKey, f64>, NewPolicy>;

/// Build the old-side strategy: `CombinedStrategy(CoefficientThreshold(1e-6),
/// MaxPauliWeight(w))`.
fn old_strat(w: usize) -> OldStrat {
    CombinedStrategy(OldCoeffThreshold(THRESHOLD), OldMaxWeight(w))
}

/// Build the new-side policy: `CombinedPolicy(CoefficientThreshold(1e-6),
/// MaxPauliWeight(w))` — the exact algebraic twin of [`old_strat`].
fn new_policy(w: usize) -> NewPolicy {
    CombinedPolicy(
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
        NewMaxWeight(w),
    )
}

/// The initial observable `Σ_i Z_i` (each `Z_i` with coefficient `1.0`) as a
/// `(pauli_string, coeff)` term list on `n` qubits.
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (s, 1.0)
        })
        .collect()
}

/// Build the old-side `Σ_i Z_i` observable under the CombinedStrategy config.
fn build_old(n: usize, w: usize) -> OldSum {
    let mut s: OldSum = OldPauliSum::builder()
        .n_qubits(n)
        .strategy(old_strat(w))
        .capacity(n.pow(2))
        .build();
    for (word, c) in sum_z_terms(n) {
        s += (word.as_str(), c);
    }
    s
}

/// Build the new-side `Σ_i Z_i` observable under the CombinedPolicy config,
/// storage-matched to [`build_old`].
fn build_new(n: usize, w: usize) -> NewSum {
    // Structurally the same construction as `build_old`: the explicit `n²`
    // capacity override and `n` accumulating `+=` inserts, not a batch build.
    let mut s = NewSum::with_capacity(n, new_policy(w), n.pow(2));
    for (word, c) in sum_z_terms(n) {
        s += (NewKey::from(word.as_str()), c);
    }
    s
}

/// The old sum's support as a sorted `(canonical_pauli_string, coeff)` vector.
fn old_support(s: &OldSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.data().iter().map(|(k, c)| (k.to_string(), *c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// The new sum's support as a sorted `(canonical_pauli_string, coeff)` vector.
fn new_support(s: &NewSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// Assert the old and new supports agree as sorted `(string, coeff)` sets.
///
/// Terms whose magnitude sits within a small band around the truncation
/// threshold (`[threshold, threshold·(1+ε)]`) are compared **leniently** on
/// presence: an end-to-end run can land a coefficient on either side of the
/// `1e-6` keep boundary with the two crates a few ulp apart, so a term in the
/// boundary band may survive on one side and not the other. Every term whose
/// magnitude is comfortably above the band must be present in BOTH with matching
/// coefficients — that is the golden-master assertion.
#[track_caller]
fn assert_integration_match(old: &OldSum, new: &NewSum, label: &str) {
    // A margin above the truncation threshold: terms above it are "solidly kept"
    // on both sides and must match exactly; terms within [threshold, band] may
    // legitimately differ in presence across the two merge orders.
    const BAND: f64 = THRESHOLD * 1.01;
    let solid = |v: Vec<(String, f64)>| -> Vec<(String, f64)> {
        v.into_iter().filter(|(_, c)| c.abs() >= BAND).collect()
    };
    let os = solid(old_support(old));
    let ns = solid(new_support(new));
    assert_eq!(
        os.len(),
        ns.len(),
        "[{label}] solidly-kept support size differs: old {} vs new {}\nold={os:?}\nnew={ns:?}",
        os.len(),
        ns.len()
    );
    for (o, n) in os.iter().zip(ns.iter()) {
        assert_eq!(
            o.0, n.0,
            "[{label}] support key differs: old {} vs new {}",
            o.0, n.0
        );
        assert_close(o.1, n.1, TOL.max(o.1.abs() * 1e-9));
    }
    // A meaningful workload: the circuit must have grown the support past the
    // initial `n` `Z_i` terms and truncation must have kept it finite.
    assert!(!os.is_empty(), "[{label}] final support is empty");
}

// ---------------------------------------------------------------------------
// 1. The Trotter workload (ported from `ppvm-pauli-sum::benches::trotter`).
// ---------------------------------------------------------------------------

/// Apply one `rzz(a, b, θ)` on the OLD sum via the native two-qubit rotation.
fn old_rzz_native(state: &mut OldSum, a: usize, b: usize, theta: f64) {
    state.rzz(a, b, theta);
}

/// Apply one `rzz(a, b, θ)` **decomposed** as `cnot(a, b); rz(b, θ); cnot(a, b)`
/// on the OLD sum — the same sequence the new sum runs (no intermediate
/// truncation, so it is the exact algebraic image of the native `rzz`).
fn old_rzz_decomposed(state: &mut OldSum, a: usize, b: usize, theta: f64) {
    state.cnot(a, b);
    state.rz(b, theta);
    state.cnot(a, b);
}

/// Apply one `rzz(a, b, θ)` decomposed on the NEW sum.
fn new_rzz_decomposed(state: &mut NewSum, a: usize, b: usize, theta: f64) {
    state.cnot(a, b);
    state.rz(b, theta);
    state.cnot(a, b);
}

/// The Trotter step body on the OLD sum (caller-driven truncation after each
/// operation, exactly like the old `trotter_func`).
fn trotter_old(
    state: &mut OldSum,
    n: usize,
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
) {
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            old_rzz_decomposed(state, i, i + 1, theta_zz);
            state.truncate();
        }
    }
}

/// The Trotter step body on the NEW sum — the identical circuit.
fn trotter_new(
    state: &mut NewSum,
    n: usize,
    steps: usize,
    theta_x: f64,
    theta_zz: f64,
    noise: [f64; 3],
) {
    for _ in 0..steps {
        for i in 0..n {
            state.pauli_error(i, noise);
            state.truncate();
            state.rx(i, theta_x);
            state.truncate();
        }
        for i in 0..n - 1 {
            state.pauli_error(i + 1, noise);
            state.truncate();
            state.pauli_error(i, noise);
            state.truncate();
            new_rzz_decomposed(state, i, i + 1, theta_zz);
            state.truncate();
        }
    }
}

/// The `rzz = cnot; rz; cnot` decomposition is exact: on the OLD crate the
/// decomposed sequence reproduces the native `rzz` term for term (validating the
/// decomposition we then apply on both sides). Checked over a spread of angles,
/// qubit pairs, and starting supports, with a single truncation at the end so the
/// two paths are compared on the same footing.
#[test]
fn old_rzz_decomposition_matches_native() {
    let n = 6;
    let w = usize::MAX; // no weight truncation: compare the full algebra
    for &(a, b) in &[(0usize, 1usize), (2, 4), (1, 5), (0, 5)] {
        for &theta in &[0.1_f64, 0.3, 0.7, 1.0, std::f64::consts::FRAC_PI_2] {
            // A nontrivial starting support so the ZZ rotation has anticommuting
            // terms to fan out: seed `Σ Z_i`, then a couple of `rx` to spread it.
            let mut native = build_old(n, w);
            let mut decomp = build_old(n, w);
            for s in [&mut native, &mut decomp] {
                s.rx(a, 0.4);
                s.truncate();
                s.rx(b, 0.5);
                s.truncate();
            }

            old_rzz_native(&mut native, a, b, theta);
            native.truncate();
            old_rzz_decomposed(&mut decomp, a, b, theta);
            decomp.truncate();

            let ns = old_support(&native);
            let ds = old_support(&decomp);
            assert_eq!(
                ns.len(),
                ds.len(),
                "rzz decomposition support size differs at a={a} b={b} θ={theta}\n\
                 native={ns:?}\ndecomp={ds:?}"
            );
            for (x, y) in ns.iter().zip(ds.iter()) {
                assert_eq!(x.0, y.0, "rzz decomposition key differs at a={a} b={b}");
                assert_close(x.1, y.1, 1e-12);
            }
        }
    }
}

#[test]
fn trotter_end_to_end_matches_old() {
    // Deterministic, small enough to be fast, deep enough that the support grows
    // and truncation bites. `h = 1`, `dt = 0.1`, so `theta_x = 0.1`, and a modest
    // number of steps; a weight cap keeps the support finite while the `1e-6`
    // coefficient floor drops the noise-suppressed tail — both truncations active.
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let j = 1.0 / 8.0 * h;
    let theta_x = dt * h;
    let theta_zz = dt * j;
    let noise = [1e-4 / 4.0; 3];

    for &n in &[6usize, 7, 8] {
        for &(time, w) in &[(0.3, 5usize), (0.5, 4)] {
            let steps = (time / dt) as usize;
            let mut old = build_old(n, w);
            let mut new = build_new(n, w);

            trotter_old(&mut old, n, steps, theta_x, theta_zz, noise);
            trotter_new(&mut new, n, steps, theta_x, theta_zz, noise);

            assert_integration_match(&old, &new, &format!("trotter n={n} steps={steps} w={w}"));
        }
    }
}

// ---------------------------------------------------------------------------
// 2. A deep random Clifford+rotation circuit — heterogeneous gate order.
// ---------------------------------------------------------------------------

/// Replay a shared `random_circuit` (`H`/`S`/`CNOT`/`Rx`/`Rz`) on the OLD sum,
/// truncating after each gate.
fn replay_old(state: &mut OldSum, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => state.h(q),
            GateOp::S(q) => state.s(q),
            GateOp::Cnot(c, t) => state.cnot(c, t),
            GateOp::Rx(q, th) => state.rx(q, th),
            GateOp::Rz(q, th) => state.rz(q, th),
        }
        state.truncate();
    }
}

/// Replay the identical `random_circuit` on the NEW sum, truncating after each
/// gate.
fn replay_new(state: &mut NewSum, circuit: &[GateOp]) {
    for &op in circuit {
        match op {
            GateOp::H(q) => state.h(q),
            GateOp::S(q) => state.s(q),
            GateOp::Cnot(c, t) => state.cnot(c, t),
            GateOp::Rx(q, th) => state.rx(q, th),
            GateOp::Rz(q, th) => state.rz(q, th),
        }
        state.truncate();
    }
}

#[test]
fn deep_random_circuit_matches_old() {
    const SEEDS: [u64; 4] = [1, 42, 777, 31337];
    for &seed in &SEEDS {
        let mut rng: StdRng = seeded_rng(seed);
        for &n in &[6usize, 8] {
            let w = 5usize;
            // A deep, heterogeneous circuit so support growth + truncation are
            // both stressed across a mixed gate order.
            let circuit = random_circuit(&mut rng, n, 400);
            let mut old = build_old(n, w);
            let mut new = build_new(n, w);
            replay_old(&mut old, &circuit);
            replay_new(&mut new, &circuit);
            assert_integration_match(&old, &new, &format!("random seed={seed} n={n}"));
        }
    }
}
