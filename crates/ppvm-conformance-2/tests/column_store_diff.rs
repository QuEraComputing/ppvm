// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! **`ColumnStore` is a drop-in alternative to `HashMapStore`** — the Phase-6
//! acceptance gate, run as a *three-way* differential.
//!
//! Every existing `-2` differential suite pins ONE new engine
//! (`Sum<HashMapStore<_, _>, _>`) against the old `ppvm-pauli-sum::PauliSum`.
//! This file replays those same workloads a **third** time on
//! `Sum<ColumnStore<_, _>, _>` — the same `Sum` engine, the same gates, the same
//! policies, the same caller-driven `truncate()` placement — and asserts the
//! observable support matches BOTH the old crate and the `HashMapStore`-backed
//! new engine. A storage backend is a representation change, so the bar is
//! observational identity, not "close enough".
//!
//! The workloads are the existing suites' verbatim, not new circuits:
//!
//! | workload | source suite |
//! |---|---|
//! | noisy-TFIM Trotter, exact + `CoefficientThreshold(1e-6)`, golden master `2.1610566562692544` | `pauli_sum_workload_diff::trotter_golden_master_and_truncation_fidelity_hold_on_both_engines` |
//! | Trotter `n = 12`, `CombinedPolicy(1e-6, MAX)`, capacity `n²` | `pauli_sum_workload_diff::trotter_n12_support_matches_old_exactly` |
//! | Trotter qubit sweep under the coefficient floor alone | `pauli_sum_workload_diff::trotter_qubit_sweep_support_and_expectation_match_old` |
//! | untruncated deep `rz;ry;rz` + CNOT-ring fan-out, plus the ℓ² norm | `pauli_sum_workload_diff::untruncated_deep_circuit_matches_old_state_and_norm` |
//! | seeded random gate stream, untruncated | `pauli_sum_workload_diff::untruncated_random_gate_stream_matches_old` |
//! | seeded 400-gate random circuit, truncated after every gate | `pauli_sum_integration_diff::deep_random_circuit_matches_old` |
//! | deferred vs eager truncation; Clifford re-key does not truncate | `pauli_sum_truncation_behaviour_diff` |
//! | exact zeros survive every gate | `pauli_sum_zero_behaviour_diff` |
//! | `preserve` snapshot-and-restore | `pauli_sum_preserve_diff` |
//! | truncation cost-grid keep-sets (weight / threshold / combined / sentinel) | `pauli_sum_truncation_boundary_diff` |
//! | pairings + zero-state contraction (GHZ backward) | `pauli_sum_pairing_order_diff`, old `tests/ghz.rs` |
//!
//! Config parity is total: `[u8; 8]` storage (`[u8; 16]` for the 128-qubit
//! truncation grid, matching that suite) and `f64` coefficients on all three
//! engines, the same policy pair, and `rzz` decomposed as
//! `cnot; rz; cnot` identically everywhere.
//!
//! **Comparison discipline.** Key sets are compared for exact equality
//! everywhere. Coefficients are compared bit-exactly where every engine applies
//! the same float operations to a term in the same order (the behaviour-contract
//! tests: zeros, deferred truncation, preserve, the truncation grid), and at a
//! tight tolerance on the deep circuits, where colliding branch keys accumulate
//! in *support order* — which is bucket order for the hash map and insertion
//! order for the columns, so the float sum legitimately reassociates. The one
//! place presence itself is compared leniently is the `1e-6` truncation band of
//! the 400-gate circuit, exactly as `pauli_sum_integration_diff` does it.

use std::collections::{BTreeMap, BTreeSet};

use ppvm_conformance_2::{GateOp, random_circuit, random_terms, seeded_rng};

// --- Old crate. ---------------------------------------------------------------
use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_word::pattern::PauliPattern as OldPattern;
use ppvm_pauli_word::word::PauliWord as OldWordT;
use ppvm_traits::traits::{
    Clifford as OldClifford, PauliError as OldPauliError, RotationOne as OldRotationOne,
    Trace as OldTrace,
};

// --- New crate: the two storage backends behind the same engine. ---------------
use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, ColumnStore, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, NoPolicy, PauliPattern as NewPattern,
    PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
    Trace as NewTrace,
};

use rand::rngs::StdRng;

/// The shared truncation floor of the integration workloads.
const THRESHOLD: f64 = 1e-6;

// --- Matched configs. `[u8; 8]` storage / `f64` coefficients on all three. -----
type NewKey = NewPauliWord<[u8; 8]>;
type OldKey = OldWordT<[u8; 8]>;
type HashSum<P> = Sum<HashMapStore<NewKey, f64>, P>;
type ColSum<P> = Sum<ColumnStore<NewKey, f64>, P>;

// --- The 128-qubit truncation grid is storage-matched to `[u8; 16]`. -----------
type NewKey16 = NewPauliWord<[u8; 16]>;
type HashSum16<P> = Sum<HashMapStore<NewKey16, f64>, P>;
type ColSum16<P> = Sum<ColumnStore<NewKey16, f64>, P>;

trait ApplyPauliNoise {
    fn apply_pauli_noise(&mut self, qubit: usize, p: [f64; 3]);
}

impl<T> ApplyPauliNoise for OldPauliSum<T>
where
    T: ppvm_traits::config::Config<Coeff = f64>,
{
    fn apply_pauli_noise(&mut self, qubit: usize, p: [f64; 3]) {
        self.pauli_error(qubit, p);
    }
}

impl<S, P> ApplyPauliNoise for Sum<S, P>
where
    S: ppvm_traits_2::Accumulate,
    S::Key: ppvm_traits_2::Word + ppvm_traits_2::Indexable,
    P: ppvm_pauli_sum_2::Policy<S::Key, S::Coeff>,
    Sum<S, P>: NewPauliError<f64>,
{
    fn apply_pauli_noise(&mut self, qubit: usize, p: [f64; 3]) {
        self.pauli_error(qubit, p, &mut ppvm_conformance_2::analytic_rng());
    }
}

// ===========================================================================
// The three-way harness.
// ===========================================================================

/// Sorted `(pauli_string, coeff)` view of an OLD sum's support.
macro_rules! old_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), *c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

/// Sorted `(pauli_string, coeff)` view of a NEW sum's support (either backend).
macro_rules! new_support {
    ($sum:expr) => {{
        let mut v: Vec<(String, f64)> = $sum.iter().map(|(k, c)| (k.to_string(), c)).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }};
}

/// The three supports produced by running ONE workload body on the old engine,
/// the `HashMapStore`-backed new engine and the `ColumnStore`-backed new engine.
struct ThreeWay {
    old: Vec<(String, f64)>,
    hash: Vec<(String, f64)>,
    column: Vec<(String, f64)>,
}

/// Run one workload body on all three engines.
///
/// The body is written once and expanded three times: the old and new crates
/// expose the same gate method names through their respective traits, and the
/// two new backends differ only in the `Sum`'s storage parameter, so a single
/// token-tree serves all three. That is what makes this a genuine *drop-in*
/// test — the driver code is literally identical.
macro_rules! three_way {
    (
        old = $old:expr,
        hash = $hash:expr,
        column = $column:expr,
        seed = $seed:expr,
        |$s:ident| $body:block
    ) => {{
        let seed: &[(String, f64)] = &$seed;
        let old = {
            let mut $s = $old;
            for (w, c) in seed {
                $s += (w.as_str(), *c);
            }
            $body
            old_support!($s)
        };
        let hash = {
            let mut $s = $hash;
            for (w, c) in seed {
                $s += (NewKey::from(w.as_str()), *c);
            }
            $body
            new_support!($s)
        };
        let column = {
            let mut $s = $column;
            for (w, c) in seed {
                $s += (NewKey::from(w.as_str()), *c);
            }
            $body
            new_support!($s)
        };
        ThreeWay { old, hash, column }
    }};
}

/// Assert two supports have the *same key set* (reported by difference) and
/// coefficients within `tol` (absolute, or relative for large coefficients).
#[track_caller]
fn assert_pair(a: &[(String, f64)], b: &[(String, f64)], tol: f64, label: &str) {
    let am: BTreeMap<&str, f64> = a.iter().map(|(k, c)| (k.as_str(), *c)).collect();
    let bm: BTreeMap<&str, f64> = b.iter().map(|(k, c)| (k.as_str(), *c)).collect();
    let only_a: Vec<&&str> = am.keys().filter(|k| !bm.contains_key(*k)).collect();
    let only_b: Vec<&&str> = bm.keys().filter(|k| !am.contains_key(*k)).collect();
    assert!(
        only_a.is_empty() && only_b.is_empty(),
        "[{label}] key sets differ: {} only left (e.g. {:?}), {} only right (e.g. {:?})",
        only_a.len(),
        only_a.first().map(|k| (k, am[**k])),
        only_b.len(),
        only_b.first().map(|k| (k, bm[**k])),
    );
    for (k, x) in &am {
        let y = bm[k];
        assert!(
            (x - y).abs() <= tol.max(x.abs() * 1e-12),
            "[{label}] coefficient at {k} differs: {x} vs {y} (tol {tol})"
        );
    }
}

/// Assert two supports are **bit-identical**: same key set, same `f64` bits
/// (so a `0.0` term must be present, and be `0.0`, in both).
#[track_caller]
fn assert_pair_exact(a: &[(String, f64)], b: &[(String, f64)], label: &str) {
    assert_eq!(
        a.len(),
        b.len(),
        "[{label}] support size differs: {} vs {}\nleft={a:?}\nright={b:?}",
        a.len(),
        b.len()
    );
    for ((ka, ca), (kb, cb)) in a.iter().zip(b.iter()) {
        assert_eq!(ka, kb, "[{label}] key differs: {ka} vs {kb}");
        assert_eq!(
            ca.to_bits(),
            cb.to_bits(),
            "[{label}] coefficient at {ka} differs: {ca} vs {cb}"
        );
    }
}

impl ThreeWay {
    /// The columnar backend matches BOTH the old crate and the hash-map backend
    /// (and, transitively, old matches hash — the standing bar, re-asserted here
    /// so a failure says which pair broke).
    #[track_caller]
    fn assert_all_match(&self, tol: f64, label: &str) {
        assert_pair(&self.old, &self.hash, tol, &format!("{label} old-vs-hash"));
        assert_pair(
            &self.old,
            &self.column,
            tol,
            &format!("{label} old-vs-column"),
        );
        assert_pair(
            &self.hash,
            &self.column,
            tol,
            &format!("{label} hash-vs-column"),
        );
        assert!(!self.old.is_empty(), "[{label}] the support is empty");
    }

    /// The stricter form: the columnar support is bit-identical to the other two.
    /// Used wherever every engine applies the same float operations to a term in
    /// the same order — in particular every zero-preservation contract.
    #[track_caller]
    fn assert_all_exact(&self, label: &str) {
        assert_pair_exact(&self.old, &self.hash, &format!("{label} old-vs-hash"));
        assert_pair_exact(&self.old, &self.column, &format!("{label} old-vs-column"));
        assert_pair_exact(&self.hash, &self.column, &format!("{label} hash-vs-column"));
    }

    /// `⟨0…0| O |0…0⟩` — the sum of the coefficients of every all-`I`/`Z` key
    /// (old's `expect_on_zero!`), per engine.
    fn zero_expectations(&self) -> (f64, f64, f64) {
        let e = |v: &Vec<(String, f64)>| -> f64 {
            v.iter()
                .filter(|(k, _)| k.chars().all(|c| c == 'I' || c == 'Z'))
                .map(|(_, c)| *c)
                .sum()
        };
        (e(&self.old), e(&self.hash), e(&self.column))
    }
}

/// `Σ_i Z_i` on `n` qubits as `(pauli_string, coeff)` terms.
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (s, 1.0)
        })
        .collect()
}

/// The noisy-TFIM first-order Trotter body with caller-driven truncation after
/// every single operation — `pauli_sum_workload_diff::trotter_evolve!` verbatim,
/// `rzz` decomposed as `cnot; rz; cnot` on every engine.
macro_rules! trotter_evolve {
    ($state:expr, $n:expr, $steps:expr, $theta_x:expr, $theta_zz:expr, $noise:expr) => {{
        for _ in 0..$steps {
            for i in 0..$n {
                $state.apply_pauli_noise(i, $noise);
                $state.truncate();
                $state.rx(i, $theta_x);
                $state.truncate();
            }
            for i in 0..$n - 1 {
                $state.apply_pauli_noise(i + 1, $noise);
                $state.truncate();
                $state.apply_pauli_noise(i, $noise);
                $state.truncate();
                $state.cnot(i, i + 1);
                $state.rz(i + 1, $theta_zz);
                $state.cnot(i, i + 1);
                $state.truncate();
            }
        }
    }};
}

// ===========================================================================
// Workload 1 — trotter-tfim: the golden master and truncation fidelity.
// ===========================================================================

/// Old's `tests/trotter.rs` acceptance bar on the COLUMNAR backend: the frozen
/// scalar `2.1610566562692544`, and the `1e-6`-floor run staying within `1e-6`
/// of the exact one (an engine that pruned at insertion drifts ~2.2e-6).
#[test]
fn trotter_golden_master_and_fidelity_hold_on_the_columnar_backend() {
    const N: usize = 4;
    const STEPS: usize = 10;
    const THETA_X: f64 = 0.1; // dt·h, dt = 0.1, h = 1
    const THETA_ZZ: f64 = 0.0125; // dt·J, J = 1/8
    const NOISE: [f64; 3] = [2.5e-5; 3]; // 1e-4 / 4 per channel
    /// `ppvm-pauli-sum/tests/trotter.rs:104`.
    const GOLDEN: f64 = 2.1610566562692544;

    let exact = three_way!(
        old = OldPauliSum::<OldByteF64<8>>::builder().n_qubits(N).build(),
        hash = HashSum::<NoPolicy>::with_policy(N, NoPolicy),
        column = ColSum::<NoPolicy>::with_policy(N, NoPolicy),
        seed = sum_z_terms(N),
        |s| {
            trotter_evolve!(s, N, STEPS, THETA_X, THETA_ZZ, NOISE);
        }
    );
    exact.assert_all_match(1e-12, "trotter N=4 exact");

    let (oe, he, ce) = exact.zero_expectations();
    for (name, v) in [("old", oe), ("hash", he), ("column", ce)] {
        assert!(
            (v - GOLDEN).abs() < 1e-9,
            "{name} exact Trotter expectation {v} drifted from golden {GOLDEN}"
        );
    }

    let approx = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(N)
            .strategy(OldCoeffThreshold(THRESHOLD))
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(
            N,
            NewCoeffThreshold {
                threshold: THRESHOLD
            }
        ),
        column = ColSum::<NewCoeffThreshold>::with_policy(
            N,
            NewCoeffThreshold {
                threshold: THRESHOLD
            }
        ),
        seed = sum_z_terms(N),
        |s| {
            trotter_evolve!(s, N, STEPS, THETA_X, THETA_ZZ, NOISE);
        }
    );
    approx.assert_all_match(1e-12, "trotter N=4 truncated");

    let (oa, ha, ca) = approx.zero_expectations();
    let drifts = [(oe - oa).abs(), (he - ha).abs(), (ce - ca).abs()];
    for (name, d) in ["old", "hash", "column"].iter().zip(drifts) {
        assert!(
            d < 1e-6,
            "{name} truncated result drifted from exact by {d:e} — an engine that \
             prunes at insertion drifts ~2.2e-6"
        );
    }
    // Sharper than "each under 1e-6": all three engines truncate identically, so
    // their drifts coincide to round-off.
    assert!(
        (drifts[0] - drifts[2]).abs() < 1e-12 && (drifts[1] - drifts[2]).abs() < 1e-12,
        "truncation drift differs across engines: {drifts:?}"
    );
}

/// The headline width: `n = 12`, 10 steps, `Combined(1e-6, MaxPauliWeight(MAX))`,
/// capacity `n²`. Strict whole-support diff plus the zero-state expectation.
#[test]
fn trotter_n12_support_matches_old_and_the_hash_backend() {
    let n = 12usize;
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let steps = ((1.0 / h) / dt) as usize;
    let theta_x = dt * h;
    let theta_zz = dt * (1.0 / 8.0 * h);
    let noise = [1e-4 / 4.0; 3];

    let old_strat = CombinedStrategy(OldCoeffThreshold(THRESHOLD), OldMaxWeight(usize::MAX));
    let new_policy = CombinedPolicy(
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
        NewMaxWeight(usize::MAX),
    );

    let out = three_way!(
        old = OldPauliSum::<OldByteF64<8, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>>::builder()
            .n_qubits(n)
            .strategy(old_strat)
            .capacity(n.pow(2))
            .build(),
        hash = HashSum::<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>::with_capacity(
            n,
            new_policy,
            n.pow(2)
        ),
        column = ColSum::<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>::with_capacity(
            n,
            new_policy,
            n.pow(2)
        ),
        seed = sum_z_terms(n),
        |s| {
            trotter_evolve!(s, n, steps, theta_x, theta_zz, noise);
        }
    );
    // Three orders of magnitude below the `1e-6` floor, so no term can be
    // reclassified kept/dropped by the tolerance.
    out.assert_all_match(1e-9, "trotter n=12");
    assert!(
        out.old.len() > 100,
        "the n=12 workload should grow a real support, got {}",
        out.old.len()
    );

    let (ov, hv, cv) = out.zero_expectations();
    let bar = 1e-9 * ov.abs().max(1.0);
    assert!(
        (ov - hv).abs() <= bar && (ov - cv).abs() <= bar,
        "zero-state expectation differs: old {ov}, hash {hv}, column {cv}"
    );
}

/// The qubit-scaling sweep under the coefficient floor ALONE. The bar is *exact*
/// support-size equality at every `n`: with a coefficient floor the surviving
/// count is a sensitive function of every accumulated coefficient, so a backend
/// that merged even slightly differently lands on a different integer.
#[test]
fn trotter_qubit_sweep_matches_on_the_columnar_backend() {
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let steps = 5usize;
    let theta_x = dt * h;
    let theta_zz = dt * 1.0; // J = 1.0: drive the support large
    let noise = [1e-4 / 4.0; 3];

    for n in (2..19usize).step_by(4) {
        let out = three_way!(
            old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
                .n_qubits(n)
                .strategy(OldCoeffThreshold(THRESHOLD))
                .capacity(n.pow(2))
                .build(),
            hash = HashSum::<NewCoeffThreshold>::with_capacity(
                n,
                NewCoeffThreshold {
                    threshold: THRESHOLD
                },
                n.pow(2)
            ),
            column = ColSum::<NewCoeffThreshold>::with_capacity(
                n,
                NewCoeffThreshold {
                    threshold: THRESHOLD
                },
                n.pow(2)
            ),
            seed = sum_z_terms(n),
            |s| {
                trotter_evolve!(s, n, steps, theta_x, theta_zz, noise);
            }
        );
        assert_eq!(
            out.old.len(),
            out.column.len(),
            "[sweep n={n}] final |support| differs: old {} vs column {}",
            out.old.len(),
            out.column.len()
        );
        assert_eq!(
            out.hash.len(),
            out.column.len(),
            "[sweep n={n}] final |support| differs: hash {} vs column {}",
            out.hash.len(),
            out.column.len()
        );
        out.assert_all_match(1e-9, &format!("sweep n={n}"));

        let (ov, hv, cv) = out.zero_expectations();
        let bar = 1e-9 * ov.abs().max(1.0);
        assert!(
            (ov - hv).abs() <= bar && (ov - cv).abs() <= bar,
            "[sweep n={n}] expectation differs: old {ov}, hash {hv}, column {cv}"
        );
    }
}

// ===========================================================================
// Workload 2 — untruncated deep circuit (pure fan-out growth).
// ===========================================================================

/// One `depth`-layer `rz(1.1); ry(2.1); rz(1.1)` + CNOT-ring circuit closed with
/// a final rotation layer (old's `benches/random-circuit.rs` shape). No
/// truncation anywhere: the support grows monotonically by fan-out, which is the
/// regime that stresses the columns' growth path (append + re-index) rather than
/// the truncation-bounded steady state.
macro_rules! random_circuit_evolve {
    ($state:expr, $n:expr, $depth:expr) => {{
        for _ in 0..$depth {
            for i in 0..$n {
                $state.rz(i, 1.1);
                $state.ry(i, 2.1);
                $state.rz(i, 1.1);
            }
            for i in 0..$n {
                $state.cnot(i, (i + 1) % $n);
            }
        }
        for i in 0..$n {
            $state.rz(i, 1.1);
            $state.ry(i, 2.1);
            $state.rz(i, 1.1);
        }
    }};
}

#[test]
fn untruncated_deep_circuit_matches_on_the_columnar_backend() {
    let n = 8usize;
    let depth = 2usize;
    let zz: String = (0..n).map(|i| if i < 2 { 'Z' } else { 'I' }).collect();

    let out = three_way!(
        old = OldPauliSum::<OldByteF64<8>>::builder().n_qubits(n).build(),
        hash = HashSum::<NoPolicy>::with_policy(n, NoPolicy),
        column = ColSum::<NoPolicy>::with_policy(n, NoPolicy),
        seed = [(zz.clone(), 1.0)],
        |s| {
            random_circuit_evolve!(s, n, depth);
        }
    );
    out.assert_all_match(1e-9, "untruncated deep circuit");
    assert!(
        out.old.len() > 500,
        "untruncated fan-out should reach a large support, got {}",
        out.old.len()
    );

    // The ℓ² norm is invariant under key ordering, so it catches a dropped or
    // duplicated term a per-key loop could miss. Computed through each new
    // engine's own `overlap` (the `Pair` impl under test), against old's support.
    let old_norm: f64 = out.old.iter().map(|(_, c)| c * c).sum();
    let mut hash_sum: HashSum<NoPolicy> = HashSum::with_policy(n, NoPolicy);
    let mut col_sum: ColSum<NoPolicy> = ColSum::with_policy(n, NoPolicy);
    hash_sum += (NewKey::from(zz.as_str()), 1.0);
    col_sum += (NewKey::from(zz.as_str()), 1.0);
    random_circuit_evolve!(hash_sum, n, depth);
    random_circuit_evolve!(col_sum, n, depth);
    let hn = hash_sum.overlap(&hash_sum);
    let cn = col_sum.overlap(&col_sum);
    let bar = 1e-9 * old_norm.abs().max(1.0);
    assert!(
        (old_norm - hn).abs() <= bar && (old_norm - cn).abs() <= bar,
        "ℓ² norm differs: old {old_norm}, hash {hn}, column {cn}"
    );
}

/// The seeded random *gate stream* (the conformance generator) with no
/// truncation: a heterogeneous `H`/`S`/`CNOT`/`Rx`/`Rz` order rather than the
/// structured layer/ring shape.
#[test]
fn untruncated_random_gate_stream_matches_on_the_columnar_backend() {
    for &seed in &[1u64, 42, 777] {
        let mut rng: StdRng = seeded_rng(seed);
        let n = 6usize;
        let circuit = random_circuit(&mut rng, n, 40);
        let zz: String = (0..n).map(|i| if i < 2 { 'Z' } else { 'I' }).collect();

        let out = three_way!(
            old = OldPauliSum::<OldByteF64<8>>::builder().n_qubits(n).build(),
            hash = HashSum::<NoPolicy>::with_policy(n, NoPolicy),
            column = ColSum::<NoPolicy>::with_policy(n, NoPolicy),
            seed = [(zz.clone(), 1.0)],
            |s| {
                for &op in &circuit {
                    match op {
                        GateOp::H(q) => s.h(q),
                        GateOp::S(q) => s.s(q),
                        GateOp::Cnot(c, t) => s.cnot(c, t),
                        GateOp::Rx(q, th) => s.rx(q, th),
                        GateOp::Rz(q, th) => s.rz(q, th),
                    }
                }
            }
        );
        out.assert_all_match(1e-9, &format!("untruncated stream seed={seed}"));
    }
}

/// The 400-gate replay with `truncate()` after EVERY gate under
/// `Combined(1e-6, MaxPauliWeight(5))` — `pauli_sum_integration_diff`'s
/// workload. Presence inside the `[τ, 1.01τ]` band is compared leniently for the
/// same reason that suite does: a term can land on either side of the keep
/// boundary with two merge orders a few ulp apart.
#[test]
fn truncated_random_circuit_replay_matches_on_the_columnar_backend() {
    /// A margin above the floor: terms above it are solidly kept everywhere.
    const BAND: f64 = THRESHOLD * 1.01;
    let solid = |v: &[(String, f64)]| -> Vec<(String, f64)> {
        v.iter().filter(|(_, c)| c.abs() >= BAND).cloned().collect()
    };

    for &seed in &[1u64, 42, 777, 31337] {
        let mut rng: StdRng = seeded_rng(seed);
        for &n in &[6usize, 8] {
            let w = 5usize;
            let circuit = random_circuit(&mut rng, n, 400);
            let old_strat = CombinedStrategy(OldCoeffThreshold(THRESHOLD), OldMaxWeight(w));
            let new_policy = CombinedPolicy(
                NewCoeffThreshold {
                    threshold: THRESHOLD,
                },
                NewMaxWeight(w),
            );

            let out = three_way!(
                    old = OldPauliSum::<
                        OldByteF64<8, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>,
                    >::builder()
                    .n_qubits(n)
                    .strategy(old_strat)
                    .capacity(n.pow(2))
                    .build(),
                    hash =
                        HashSum::<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>::with_capacity(
                            n,
                            new_policy,
                            n.pow(2)
                        ),
                    column =
                        ColSum::<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>::with_capacity(
                            n,
                            new_policy,
                            n.pow(2)
                        ),
                    seed = sum_z_terms(n),
                    |s| {
                        for &op in &circuit {
                            match op {
                                GateOp::H(q) => s.h(q),
                                GateOp::S(q) => s.s(q),
                                GateOp::Cnot(c, t) => s.cnot(c, t),
                                GateOp::Rx(q, th) => s.rx(q, th),
                                GateOp::Rz(q, th) => s.rz(q, th),
                            }
                            s.truncate();
                        }
                    }
                );
            let label = format!("truncated replay seed={seed} n={n}");
            let (o, h, c) = (solid(&out.old), solid(&out.hash), solid(&out.column));
            assert!(!o.is_empty(), "[{label}] final support is empty");
            assert_pair(&o, &h, 1e-9, &format!("{label} old-vs-hash"));
            assert_pair(&o, &c, 1e-9, &format!("{label} old-vs-column"));
            assert_pair(&h, &c, 1e-9, &format!("{label} hash-vs-column"));
        }
    }
}

// ===========================================================================
// Behaviour contracts — WHEN side effects fire, on the columnar backend.
// ===========================================================================

/// Contract 1 (truncation timing): gates NEVER truncate. Two `rx(0.03)` with no
/// truncate in between merge into a `Y ≈ sin(2θ) ≈ 0.06` above the `τ = 0.05`
/// floor, which the final truncate must keep — on all three engines.
#[test]
fn deferred_truncation_matches_on_the_columnar_backend() {
    const TAU: f64 = 0.05;
    const THETA: f64 = 0.03;

    let deferred = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(1)
            .strategy(OldCoeffThreshold(TAU))
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(1, NewCoeffThreshold { threshold: TAU }),
        column = ColSum::<NewCoeffThreshold>::with_policy(1, NewCoeffThreshold { threshold: TAU }),
        seed = [("Z".to_string(), 1.0)],
        |s| {
            s.rx(0, THETA);
            s.rx(0, THETA);
            s.truncate();
        }
    );
    assert!(
        deferred.old.iter().any(|(k, c)| k == "Y" && c.abs() >= TAU),
        "test setup broken: old should keep an above-threshold Y, got {:?}",
        deferred.old
    );
    deferred.assert_all_exact("two rx, truncate deferred");

    // The eager schedule genuinely loses the sub-threshold Y — every engine must
    // agree on THAT answer too (the fix removed the internal truncate without
    // disturbing the explicit one).
    let eager = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(1)
            .strategy(OldCoeffThreshold(TAU))
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(1, NewCoeffThreshold { threshold: TAU }),
        column = ColSum::<NewCoeffThreshold>::with_policy(1, NewCoeffThreshold { threshold: TAU }),
        seed = [("Z".to_string(), 1.0)],
        |s| {
            s.rx(0, THETA);
            s.truncate();
            s.rx(0, THETA);
            s.truncate();
        }
    );
    assert!(
        !eager.old.iter().any(|(k, _)| k == "Y"),
        "test setup broken: eager truncation should have dropped the Y, got {:?}",
        eager.old
    );
    eager.assert_all_exact("two rx, truncate after each");
}

/// Contract 3 (re-key does not truncate): under `MaxPauliWeight(1)` a
/// `cnot(0, 1)` maps `XI ↦ XX` (weight 2) and the over-weight term must STAY
/// until an explicit `truncate()` — including on the columnar re-key, which
/// rebuilds its planes in place.
#[test]
fn clifford_rekey_does_not_truncate_on_the_columnar_backend() {
    let after_gate = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldMaxWeight>>::builder()
            .n_qubits(2)
            .strategy(OldMaxWeight(1))
            .build(),
        hash = HashSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        column = ColSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        seed = [("XI".to_string(), 1.0)],
        |s| {
            s.cnot(0, 1);
        }
    );
    assert!(
        after_gate.old.iter().any(|(k, _)| k == "XX"),
        "test setup broken: old should still hold the over-weight XX"
    );
    after_gate.assert_all_exact("after cnot, no truncate");

    // The gate is an involution: a round trip is the identity on every engine.
    let round_trip = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldMaxWeight>>::builder()
            .n_qubits(2)
            .strategy(OldMaxWeight(1))
            .build(),
        hash = HashSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        column = ColSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        seed = [("XI".to_string(), 1.0)],
        |s| {
            s.cnot(0, 1);
            s.cnot(0, 1);
        }
    );
    assert_eq!(round_trip.column, vec![("XI".to_string(), 1.0)]);
    round_trip.assert_all_exact("cnot round trip");

    // …and the explicit truncate still fires.
    let truncated = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldMaxWeight>>::builder()
            .n_qubits(2)
            .strategy(OldMaxWeight(1))
            .build(),
        hash = HashSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        column = ColSum::<NewMaxWeight>::with_policy(2, NewMaxWeight(1)),
        seed = [("XI".to_string(), 1.0)],
        |s| {
            s.cnot(0, 1);
            s.truncate();
        }
    );
    assert!(
        truncated.column.is_empty(),
        "explicit truncate should have dropped the over-weight XX, got {:?}",
        truncated.column
    );
    truncated.assert_all_exact("cnot then explicit truncate");
}

/// Contract 2 (no implicit reduce): an exactly-zero coefficient stays a live
/// entry. This is the sharpest available test of the columnar backend, whose
/// `reduce` IS a prefix-sum compaction — it must run only when the caller asks.
#[test]
fn exact_zeros_survive_every_gate_on_the_columnar_backend() {
    const SEEDS: [u64; 4] = [1, 42, 777, 31337];
    const WIDTHS: [usize; 4] = [1, 3, 5, 8];

    // (a) The identity rotation adds a 0.0-coefficient branch key.
    let ident = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(1)
            .strategy(OldCoeffThreshold(1e-12))
            .build(),
        hash = HashSum::<NoPolicy>::with_policy(1, NoPolicy),
        column = ColSum::<NoPolicy>::with_policy(1, NoPolicy),
        seed = [("Z".to_string(), 1.0)],
        |s| {
            s.rx(0, 0.0);
        }
    );
    assert_eq!(ident.old.len(), 2, "old keeps the 0.0 Y branch");
    ident.assert_all_exact("identity rotation zero branch");

    for &seed in &SEEDS {
        let mut rng = seeded_rng(seed);
        for &w in &WIDTHS {
            let terms = random_terms(&mut rng, w, 12);
            for q in 0..w {
                // (b) A zero channel eigenvalue: λ_X = 1 − 2(p_Y + p_Z) = 0.
                let ch = three_way!(
                    old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
                        .n_qubits(w)
                        .strategy(OldCoeffThreshold(1e-12))
                        .build(),
                    hash = HashSum::<NoPolicy>::with_policy(w, NoPolicy),
                    column = ColSum::<NoPolicy>::with_policy(w, NoPolicy),
                    seed = terms.clone(),
                    |s| {
                        s.apply_pauli_noise(q, [0.0, 0.25, 0.25]);
                    }
                );
                ch.assert_all_exact(&format!("zero eigenvalue seed={seed} w={w} q={q}"));

                // (c) The identity rotation on a random support.
                let rot = three_way!(
                    old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
                        .n_qubits(w)
                        .strategy(OldCoeffThreshold(1e-12))
                        .build(),
                    hash = HashSum::<NoPolicy>::with_policy(w, NoPolicy),
                    column = ColSum::<NoPolicy>::with_policy(w, NoPolicy),
                    seed = terms.clone(),
                    |s| {
                        s.ry(q, 0.0);
                    }
                );
                rot.assert_all_exact(&format!("zero branch seed={seed} w={w} q={q}"));
            }

            // (d) `*= 0.0` keeps the whole key set (old's `test_reset_channel`
            //     shape) — the columnar `scale` is one contiguous pass and must
            //     not compact.
            let scaled = three_way!(
                old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
                    .n_qubits(w)
                    .strategy(OldCoeffThreshold(1e-12))
                    .build(),
                hash = HashSum::<NoPolicy>::with_policy(w, NoPolicy),
                column = ColSum::<NoPolicy>::with_policy(w, NoPolicy),
                seed = terms.clone(),
                |s| {
                    s *= 0.0;
                }
            );
            assert!(
                scaled.column.iter().all(|(_, c)| *c == 0.0),
                "`*= 0.0` must zero every coefficient"
            );
            scaled.assert_all_exact(&format!("scale by zero seed={seed} w={w}"));
        }
    }

    // (e) An inserted 0.0 and an exact cancellation both keep their key.
    for &w in &WIDTHS {
        let key = "Z".repeat(w);
        let inserted = three_way!(
            old = OldPauliSum::<OldByteF64<8>>::builder().n_qubits(w).build(),
            hash = HashSum::<NoPolicy>::with_policy(w, NoPolicy),
            column = ColSum::<NoPolicy>::with_policy(w, NoPolicy),
            seed = [(key.clone(), 0.0)],
            |s| {
                let _ = &s;
            }
        );
        assert_eq!(inserted.column.len(), 1, "an inserted 0.0 keeps its key");
        inserted.assert_all_exact("inserted zero");

        let cancelled = three_way!(
            old = OldPauliSum::<OldByteF64<8>>::builder().n_qubits(w).build(),
            hash = HashSum::<NoPolicy>::with_policy(w, NoPolicy),
            column = ColSum::<NoPolicy>::with_policy(w, NoPolicy),
            seed = [(key.clone(), 1.5), (key.clone(), -1.5)],
            |s| {
                let _ = &s;
            }
        );
        assert_eq!(cancelled.column, vec![(key.clone(), 0.0)]);
        cancelled.assert_all_exact("exact cancellation");
    }
}

/// Contract 5 (`preserve` snapshot-and-restore) on the columnar backend: the
/// preserved keys the policy drops come back at their PRE-truncate coefficient,
/// a preserved key the policy KEPT is not doubled, and an absent preserved key
/// is not conjured into the support.
#[test]
fn preserve_set_restore_matches_on_the_columnar_backend() {
    const TAU: f64 = 0.5;
    let single_z = |n: usize| -> Vec<String> {
        (0..n)
            .map(|i| (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect())
            .collect()
    };

    // (a) The coefficient floor drops the preserved keys; they are restored.
    let terms: Vec<(String, f64)> = [
        ("ZII", 1e-6),
        ("IZI", 1e-6),
        ("IIZ", 1e-6),
        ("XYZ", 1e-6),
        ("XXX", 0.7),
    ]
    .iter()
    .map(|(s, c)| (s.to_string(), *c))
    .collect();

    let out = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(3)
            .strategy(OldCoeffThreshold(TAU))
            .preserve_strings(single_z(3).into_iter().map(OldKey::from).collect())
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(3, NewCoeffThreshold { threshold: TAU })
            .preserving(single_z(3).into_iter().map(|s| NewKey::from(s.as_str()))),
        column = ColSum::<NewCoeffThreshold>::with_policy(3, NewCoeffThreshold { threshold: TAU })
            .preserving(single_z(3).into_iter().map(|s| NewKey::from(s.as_str()))),
        seed = terms,
        |s| {
            s.truncate();
        }
    );
    out.assert_all_exact("coefficient floor + preserve set");
    let keys: Vec<&str> = out.column.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, vec!["IIZ", "IZI", "XXX", "ZII"]);

    // (b) Repeated propagate + truncate: restored every round, at the
    //     pre-truncate coefficient.
    const N: usize = 4;
    const THETA: f64 = 0.37;
    const ROUNDS: usize = 12;
    let repeated = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(N)
            .strategy(OldCoeffThreshold(TAU))
            .preserve_strings(single_z(N).into_iter().map(OldKey::from).collect())
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(N, NewCoeffThreshold { threshold: TAU })
            .preserving(single_z(N).into_iter().map(|s| NewKey::from(s.as_str()))),
        column = ColSum::<NewCoeffThreshold>::with_policy(N, NewCoeffThreshold { threshold: TAU })
            .preserving(single_z(N).into_iter().map(|s| NewKey::from(s.as_str()))),
        seed = single_z(N)
            .into_iter()
            .map(|s| (s, 1.0))
            .collect::<Vec<_>>(),
        |s| {
            for _ in 0..ROUNDS {
                for q in 0..N {
                    s.rx(q, THETA);
                }
                s.truncate();
            }
        }
    );
    repeated.assert_all_exact("repeated rx + truncate under a preserve set");
    for s in single_z(N) {
        assert!(
            repeated.column.iter().any(|(k, _)| *k == s),
            "preserved key {s} must survive on the columnar backend"
        );
    }

    // (c) Composes with the weight cap; a kept key is never double-added; an
    //     absent preserved key is not inserted.
    let weight_terms: Vec<(String, f64)> = [("ZII", 2.0), ("ZZI", 3.0), ("XYZ", 5.0)]
        .iter()
        .map(|(s, c)| (s.to_string(), *c))
        .collect();
    let capped = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldMaxWeight>>::builder()
            .n_qubits(3)
            .strategy(OldMaxWeight(1))
            .preserve_strings(["ZII", "ZZI"].into_iter().map(OldKey::from).collect())
            .build(),
        hash = HashSum::<NewMaxWeight>::with_policy(3, NewMaxWeight(1))
            .preserving(["ZII", "ZZI"].into_iter().map(NewKey::from)),
        column = ColSum::<NewMaxWeight>::with_policy(3, NewMaxWeight(1))
            .preserving(["ZII", "ZZI"].into_iter().map(NewKey::from)),
        seed = weight_terms,
        |s| {
            s.truncate();
        }
    );
    capped.assert_all_exact("weight cap + preserve set");
    assert_eq!(
        capped.column,
        vec![("ZII".to_string(), 2.0), ("ZZI".to_string(), 3.0)],
        "a preserved key the policy KEPT must not be doubled"
    );

    let absent = three_way!(
        old = OldPauliSum::<OldByteF64<8, OldCoeffThreshold>>::builder()
            .n_qubits(2)
            .strategy(OldCoeffThreshold(TAU))
            .preserve_strings(["ZI", "IZ"].into_iter().map(OldKey::from).collect())
            .build(),
        hash = HashSum::<NewCoeffThreshold>::with_policy(2, NewCoeffThreshold { threshold: TAU })
            .preserving(["ZI", "IZ"].into_iter().map(NewKey::from)),
        column = ColSum::<NewCoeffThreshold>::with_policy(2, NewCoeffThreshold { threshold: TAU })
            .preserving(["ZI", "IZ"].into_iter().map(NewKey::from)),
        seed = [("ZI".to_string(), 1e-9)],
        |s| {
            s.truncate();
        }
    );
    absent.assert_all_exact("absent preserved key");
    assert_eq!(
        absent.column.len(),
        1,
        "an absent preserved key is not created"
    );
}

// ===========================================================================
// The truncation cost-grid keep-sets (128 qubits, `[u8; 16]` storage).
// ===========================================================================

/// Qubit count of the truncation workload (old's `benches/truncation-weight.rs`).
const GRID_N: usize = 128;
/// Term count.
const GRID_TERMS: usize = 1000;
/// The coefficient floor used by the grid's threshold cells.
const GRID_TAU: f64 = 1e-12;

/// `pauli_sum_truncation_boundary_diff::profile_terms`, verbatim.
fn profile_terms(target_weight: usize) -> Vec<(String, f64)> {
    let stride = (GRID_N / target_weight).max(1);
    (0..GRID_TERMS)
        .map(|k| {
            let mut w = vec!['I'; GRID_N];
            for j in 0..target_weight {
                let pos = (j * stride + k) % GRID_N;
                w[pos] = ['X', 'Y', 'Z'][(k + j) % 3];
            }
            let extra = (k * 7 + 3) % GRID_N;
            if w[extra] == 'I' {
                w[extra] = ['X', 'Y', 'Z'][k % 3];
            }
            (w.into_iter().collect::<String>(), 1.0 / (k as f64 + 1.0))
        })
        .collect()
}

/// A word with exactly `weight` non-identity sites (packed at the front).
fn word_of_weight(weight: usize, letter: char) -> String {
    (0..GRID_N)
        .map(|i| if i < weight { letter } else { 'I' })
        .collect()
}

/// The boundary terms: exactly at and just under the coefficient floor, and
/// exactly at and just over each weight cutoff under test.
fn boundary_terms(cutoffs: &[usize]) -> Vec<(String, f64)> {
    let mut v = vec![
        (word_of_weight(1, 'X'), GRID_TAU), // |c| == τ → KEPT
        (
            word_of_weight(2, 'X'),
            f64::from_bits(GRID_TAU.to_bits() - 1),
        ), // τ − 1ulp → dropped
        (
            word_of_weight(3, 'X'),
            f64::from_bits(GRID_TAU.to_bits() + 1),
        ), // τ + 1ulp → kept
        (word_of_weight(4, 'X'), -GRID_TAU), // the rule is on the magnitude
    ];
    for &w in cutoffs {
        if w == 0 || w >= GRID_N {
            continue;
        }
        v.push((word_of_weight(w, 'Y'), 1.0)); // weight == w → KEPT
        v.push((word_of_weight(w + 1, 'Y'), 1.0)); // weight == w + 1 → dropped
    }
    v
}

macro_rules! old_grid_survivors {
    ($cfg:ty, $strat:expr, $terms:expr) => {{
        let mut s: OldPauliSum<$cfg> = OldPauliSum::builder()
            .n_qubits(GRID_N)
            .strategy($strat)
            .capacity(GRID_TERMS * 2)
            .build();
        for (w, c) in $terms {
            s += (w.as_str(), *c);
        }
        s.truncate();
        s.data()
            .iter()
            .map(|(k, _)| k.to_string())
            .collect::<BTreeSet<String>>()
    }};
}

macro_rules! new_grid_survivors {
    ($sum_ty:ty, $policy:expr, $terms:expr) => {{
        let mut s: $sum_ty = Sum::with_capacity(GRID_N, $policy, GRID_TERMS * 2);
        for (w, c) in $terms {
            s += (NewKey16::from(w.as_str()), *c);
        }
        s.truncate();
        s.iter()
            .map(|(k, _)| k.to_string())
            .collect::<BTreeSet<String>>()
    }};
}

#[track_caller]
fn assert_same_keys(a: &BTreeSet<String>, b: &BTreeSet<String>, label: &str) {
    let only_a: Vec<&String> = a.difference(b).collect();
    let only_b: Vec<&String> = b.difference(a).collect();
    assert!(
        only_a.is_empty() && only_b.is_empty(),
        "[{label}] keep-sets differ: {} only left, {} only right",
        only_a.len(),
        only_b.len()
    );
}

/// Contracts 6 & 7 on the columnar backend: the `MaxPauliWeight` / threshold /
/// combined keep-rules and the `usize::MAX` disable sentinel, over the old
/// truncation cost grid. The columnar `Retain` is a *compaction*, a different
/// mechanism from the hash map's `retain`, so the boundary cells are the direct
/// detector for an off-by-one there.
#[test]
fn truncation_grid_keep_sets_match_on_the_columnar_backend() {
    const CUTOFFS: [usize; 4] = [10, 100, 1000, usize::MAX];
    for target in [3usize, 50, 120] {
        let mut terms = profile_terms(target);
        terms.extend(boundary_terms(&CUTOFFS));

        // --- MaxPauliWeight cells, including the disable sentinel. ----------
        for w in CUTOFFS {
            let old = old_grid_survivors!(OldByteF64<16, OldMaxWeight>, OldMaxWeight(w), &terms);
            let hash = new_grid_survivors!(HashSum16<NewMaxWeight>, NewMaxWeight(w), &terms);
            let col = new_grid_survivors!(ColSum16<NewMaxWeight>, NewMaxWeight(w), &terms);
            let label = format!("weight profile={target} cutoff={w}");
            assert_same_keys(&old, &hash, &format!("{label} old-vs-hash"));
            assert_same_keys(&old, &col, &format!("{label} old-vs-column"));
            assert_same_keys(&hash, &col, &format!("{label} hash-vs-column"));
            if w == usize::MAX {
                let all: BTreeSet<String> = terms.iter().map(|(k, _)| k.clone()).collect();
                assert_eq!(
                    col, all,
                    "the sentinel dropped a term on the columnar backend"
                );
            } else {
                for k in &col {
                    let weight = k.chars().filter(|&c| c != 'I').count();
                    assert!(weight <= w, "kept a weight-{weight} term under cutoff {w}");
                }
            }
        }

        // --- CoefficientThreshold cell: `|c| >= τ` keeps. -------------------
        let mut tterms = profile_terms(target);
        tterms.extend(boundary_terms(&[]));
        let old = old_grid_survivors!(
            OldByteF64<16, OldCoeffThreshold>,
            OldCoeffThreshold(GRID_TAU),
            &tterms
        );
        let hash = new_grid_survivors!(
            HashSum16<NewCoeffThreshold>,
            NewCoeffThreshold {
                threshold: GRID_TAU
            },
            &tterms
        );
        let col = new_grid_survivors!(
            ColSum16<NewCoeffThreshold>,
            NewCoeffThreshold {
                threshold: GRID_TAU
            },
            &tterms
        );
        let label = format!("threshold profile={target}");
        assert_same_keys(&old, &hash, &format!("{label} old-vs-hash"));
        assert_same_keys(&old, &col, &format!("{label} old-vs-column"));
        assert_same_keys(&hash, &col, &format!("{label} hash-vs-column"));
        assert!(
            col.contains(&word_of_weight(1, 'X')),
            "|c| == τ must be KEPT on the columnar backend"
        );
        assert!(
            !col.contains(&word_of_weight(2, 'X')),
            "|c| == τ − 1ulp must be dropped on the columnar backend"
        );
        assert!(
            col.contains(&word_of_weight(4, 'X')),
            "|c| == −τ must be KEPT (the rule is on the magnitude)"
        );

        // --- Combined cells: two sequential retain passes. -------------------
        for w in [10usize, 100, usize::MAX] {
            let old = old_grid_survivors!(
                OldByteF64<16, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>,
                CombinedStrategy(OldCoeffThreshold(GRID_TAU), OldMaxWeight(w)),
                &terms
            );
            let hash = new_grid_survivors!(
                HashSum16<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>,
                CombinedPolicy(
                    NewCoeffThreshold {
                        threshold: GRID_TAU
                    },
                    NewMaxWeight(w)
                ),
                &terms
            );
            let col = new_grid_survivors!(
                ColSum16<CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>,
                CombinedPolicy(
                    NewCoeffThreshold {
                        threshold: GRID_TAU
                    },
                    NewMaxWeight(w)
                ),
                &terms
            );
            let label = format!("combined profile={target} cutoff={w}");
            assert_same_keys(&old, &hash, &format!("{label} old-vs-hash"));
            assert_same_keys(&old, &col, &format!("{label} old-vs-column"));
            assert_same_keys(&hash, &col, &format!("{label} hash-vs-column"));
        }
    }
}

/// Contract 14 (`capacity()` reports the resolved hint) is backend-independent
/// by construction, but the columnar store must actually *size both buffers*
/// from it — so assert the hint survives the swap.
#[test]
fn capacity_hint_is_reported_identically_by_both_backends() {
    for n in [0usize, 1, 4, 12, 128] {
        let policy = NewCoeffThreshold { threshold: 1e-12 };
        let hash: HashSum<NewCoeffThreshold> = HashSum::with_policy(n, policy);
        let col: ColSum<NewCoeffThreshold> = ColSum::with_policy(n, policy);
        assert_eq!(hash.capacity(), col.capacity(), "policy hint at n={n}");
        assert_eq!(col.capacity(), n * 10, "the resolved hint is n*10");

        let hash: HashSum<NewCoeffThreshold> = HashSum::with_capacity(n, policy, 77);
        let col: ColSum<NewCoeffThreshold> = ColSum::with_capacity(n, policy, 77);
        assert_eq!(hash.capacity(), 77);
        assert_eq!(col.capacity(), 77, "the explicit override must be reported");
        assert_eq!(col.n_sites(), n);
    }
}

// ===========================================================================
// Pairings and the zero-state contraction.
// ===========================================================================

/// Contract 15 (trace / overlap) on the columnar backend, over the pairing suite's
/// cases plus old's frozen GHZ-backward scalar.
#[test]
fn pairings_and_zero_state_trace_match_on_the_columnar_backend() {
    /// One pairing case: the two term lists to build `a` and `b` from.
    type Case = (
        &'static [(&'static str, f64)],
        &'static [(&'static str, f64)],
    );
    // (a) The pairing cases: empty, self, orthogonal, partial overlap, zeros.
    let cases: [Case; 6] = [
        (&[("IIII", 1.0), ("XIII", 2.0)], &[]),
        (&[("XIII", 3.0)], &[("XIII", 3.0)]),
        (&[("XIII", 1.0)], &[("ZIII", 1.0)]),
        (
            &[("XIII", 1.0), ("YIII", 2.0), ("ZIII", 3.0)],
            &[("YIII", 5.0), ("ZIII", 7.0), ("IIII", 11.0)],
        ),
        (&[("XIII", 1.0), ("IIII", 0.0)], &[("IIII", 9.0)]),
        (&[], &[("IIII", 1.0)]),
    ];
    for (i, (a, b)) in cases.iter().enumerate() {
        let build_old = |t: &[(&str, f64)]| {
            let mut s: OldPauliSum<OldByteF64<8>> = OldPauliSum::builder().n_qubits(4).build();
            for (w, c) in t {
                s += (*w, *c);
            }
            s
        };
        let mut ha: HashSum<NoPolicy> = HashSum::with_policy(4, NoPolicy);
        let mut hb: HashSum<NoPolicy> = HashSum::with_policy(4, NoPolicy);
        let mut ca: ColSum<NoPolicy> = ColSum::with_policy(4, NoPolicy);
        let mut cb: ColSum<NoPolicy> = ColSum::with_policy(4, NoPolicy);
        for (w, c) in a.iter() {
            ha += (NewKey::from(*w), *c);
            ca += (NewKey::from(*w), *c);
        }
        for (w, c) in b.iter() {
            hb += (NewKey::from(*w), *c);
            cb += (NewKey::from(*w), *c);
        }
        let o = build_old(a).overlap(&build_old(b));
        let h = ha.overlap(&hb);
        let cc = ca.overlap(&cb);
        let bar = 1e-12 * o.abs().max(1.0);
        assert!(
            (o - h).abs() <= bar && (o - cc).abs() <= bar,
            "[pairing case {i}] overlap differs: old {o}, hash {h}, column {cc}"
        );
        // Symmetry, on the columnar backend too.
        assert!((cc - cb.overlap(&ca)).abs() <= bar, "overlap is symmetric");
    }

    // (b) The frozen GHZ-backward scalar (`ppvm-pauli-sum/tests/ghz.rs`):
    //     seed ZZ, cnot(0, 1), h(0), contract against the zero-state pattern.
    let mut old: OldPauliSum<OldByteF64<8>> = OldPauliSum::builder().n_qubits(2).build();
    old += ("ZZ", 1.0);
    old.cnot(0, 1);
    old.h(0);
    let old_pattern: OldPattern = "Z?*".into();
    let ot: f64 = old.trace(&old_pattern);

    let mut hash: HashSum<NoPolicy> = HashSum::with_policy(2, NoPolicy);
    let mut col: ColSum<NoPolicy> = ColSum::with_policy(2, NoPolicy);
    hash += (NewKey::from("ZZ"), 1.0);
    col += (NewKey::from("ZZ"), 1.0);
    hash.cnot(0, 1);
    hash.h(0);
    col.cnot(0, 1);
    col.h(0);
    let pattern = NewPattern::zero_state();
    let ht: f64 = hash.trace(&pattern);
    let ct: f64 = col.trace(&pattern);
    assert_eq!(ot, 1.0, "old's frozen GHZ-backward scalar");
    assert_eq!(ht, 1.0, "hash backend GHZ-backward scalar");
    assert_eq!(ct, 1.0, "columnar backend GHZ-backward scalar");
}
