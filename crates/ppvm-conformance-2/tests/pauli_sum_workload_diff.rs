// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The **integration-baseline numeric bars**: the old crate's real workloads,
//! replayed end to end on BOTH engines under an identical algebraic config, with
//! the acceptance thresholds the old crate itself froze.
//!
//! `pauli_sum_integration_diff.rs` already ports the *shape* of the Trotter and
//! random-circuit workloads and diffs the two supports leniently around the
//! truncation boundary. This file adds the parts of the acceptance bar that
//! leniency cannot express:
//!
//! 1. **Golden master** (`ppvm-pauli-sum/tests/trotter.rs`): the exact,
//!    untruncated `N = 4`, `STEPS = 10` noisy-TFIM Trotter expectation must equal
//!    the frozen constant `2.1610566562692544` to `1e-9`, on BOTH engines. A
//!    frozen constant catches a sign/addressing/ordering error that an old-vs-new
//!    diff cannot: a bug present in both engines still diffs clean.
//! 2. **Truncation fidelity**: the `CoefficientThreshold(1e-6)` run must stay
//!    within `1e-6` of the exact run. Old drifts ~2.4e-8; an engine that prunes
//!    *at insertion* (instead of deferring truncation to the caller) drifts
//!    ~2.2e-6 and fails here — this is the numeric form of the deferred-truncation
//!    contract.
//! 3. **Strict** end-to-end support diff at the headline width `n = 12`: identical
//!    key sets (exact set equality, not a size check) and per-key coefficients to
//!    `1e-9`, three orders of magnitude below the `1e-6` floor.
//! 4. **Qubit-scaling sweep**: the same circuit over a range of `n` under
//!    `CoefficientThreshold` alone, asserting the final `|support|` is *exactly*
//!    equal old-vs-new at every `n` (the coefficient floor makes support size a
//!    sensitive function of every accumulated coefficient) and the zero-state
//!    expectation agrees to `1e-9` relative.
//! 5. **Untruncated deep circuit** (`NoStrategy` / `NoPolicy`): support grows by
//!    pure fan-out, so the whole state is diffed strictly — plus `overlap(f, f)`
//!    (the ℓ² norm), which is invariant under key ordering and catches a dropped
//!    or duplicated term a per-key loop with a size check could miss.
//!
//! Config parity is total: `[u8; 8]` storage and `f64` coefficients on both
//! sides, the same policy/strategy pair, the same caller-driven `truncate()`
//! placement, and the same `rzz(a, b, θ) = cnot(a, b); rz(b, θ); cnot(a, b)`
//! decomposition (validated against old's native `rzz` in
//! `pauli_sum_integration_diff::old_rzz_decomposition_matches_native`).

use ppvm_conformance_2::seeded_rng;

use ppvm_pauli_sum::config::fxhash::ByteF64 as OldByteF64;
use ppvm_pauli_sum::strategy::{
    CoefficientThreshold as OldCoeffThreshold, CombinedStrategy, MaxPauliWeight as OldMaxWeight,
};
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_traits::traits::{
    Clifford as OldClifford, PauliError as OldPauliError, RotationOne as OldRotationOne,
};

use ppvm_pauli_sum_2::{
    CoefficientThreshold as NewCoeffThreshold, CombinedPolicy, HashMapStore,
    MaxPauliWeight as NewMaxWeight, NoPolicy, PauliWord as NewPauliWord, Sum,
};
use ppvm_traits_2::{
    Clifford as NewClifford, PauliError as NewPauliError, RotationOne as NewRotationOne,
};

use rand::rngs::StdRng;
use std::collections::BTreeMap;

/// The shared truncation floor.
const THRESHOLD: f64 = 1e-6;

// --- Matched configs. Storage `[u8; 8]` (64-qubit capacity) on both sides. ----
type NewKey = NewPauliWord<[u8; 8]>;

/// Exact (untruncated) pair.
type OldExactSum = OldPauliSum<OldByteF64<8>>;
type NewExactSum = Sum<HashMapStore<NewKey, f64>, NoPolicy>;

/// Coefficient-floor-only pair (the qubit-sweep and fidelity configs).
type OldThreshSum = OldPauliSum<OldByteF64<8, OldCoeffThreshold>>;
type NewThreshSum = Sum<HashMapStore<NewKey, f64>, NewCoeffThreshold>;

/// The headline pair: `Combined(CoefficientThreshold(1e-6), MaxPauliWeight(MAX))`.
type OldCombinedSum = OldPauliSum<OldByteF64<8, CombinedStrategy<OldCoeffThreshold, OldMaxWeight>>>;
type NewCombinedSum =
    Sum<HashMapStore<NewKey, f64>, CombinedPolicy<NewCoeffThreshold, NewMaxWeight>>;

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

/// `Σ_i Z_i` as `(pauli_string, coeff)` terms on `n` qubits.
fn sum_z_terms(n: usize) -> Vec<(String, f64)> {
    (0..n)
        .map(|i| {
            let s: String = (0..n).map(|j| if j == i { 'Z' } else { 'I' }).collect();
            (s, 1.0)
        })
        .collect()
}

/// Seed an OLD sum with `Σ_i Z_i` through the `+=` string-add path (old's own
/// construction in `tests/trotter.rs` and `benches/trotter.rs`).
macro_rules! seed_old {
    ($state:expr, $n:expr) => {{
        for (w, c) in sum_z_terms($n) {
            $state += (w.as_str(), c);
        }
    }};
}

/// Seed a NEW sum with `Σ_i Z_i` through the same accumulating `+=` path.
macro_rules! seed_new {
    ($state:expr, $n:expr) => {{
        for (w, c) in sum_z_terms($n) {
            $state += (NewKey::from(w.as_str()), c);
        }
    }};
}

/// The noisy-TFIM first-order Trotter body, caller-driven truncation after every
/// single operation — old's `trotter_func` verbatim, with `rzz` decomposed as
/// `cnot; rz; cnot` so the *same* gate sequence runs on both engines.
///
/// A macro rather than a generic fn: the two engines expose the same method
/// names through different traits, so one body serves both without replaying
/// each crate's bound list.
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

/// `⟨0…0| O |0…0⟩` on the OLD engine: the sum of the coefficients of every
/// all-`I`/`Z` key (old's `expect_on_zero!`).
macro_rules! old_zero_expectation {
    ($state:expr) => {{
        let mut acc = 0.0_f64;
        for (k, v) in $state.data().iter() {
            if k.to_string().chars().all(|c| c == 'I' || c == 'Z') {
                acc += *v;
            }
        }
        acc
    }};
}

/// `⟨0…0| O |0…0⟩` on the NEW engine — the identical contraction.
macro_rules! new_zero_expectation {
    ($state:expr) => {{
        let mut acc = 0.0_f64;
        for (k, v) in $state.iter() {
            if k.to_string().chars().all(|c| c == 'I' || c == 'Z') {
                acc += v;
            }
        }
        acc
    }};
}

/// The OLD support as a key→coeff map.
macro_rules! old_support {
    ($state:expr) => {{
        let m: BTreeMap<String, f64> = $state
            .data()
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        m
    }};
}

/// The NEW support as a key→coeff map.
macro_rules! new_support {
    ($state:expr) => {{
        let m: BTreeMap<String, f64> = $state.iter().map(|(k, v)| (k.to_string(), v)).collect();
        m
    }};
}

/// Assert two supports are **strictly** identical: same key set (exact set
/// equality, reported by difference rather than by size) and per-key
/// coefficients within `tol`.
#[track_caller]
fn assert_supports_identical(
    old: &BTreeMap<String, f64>,
    new: &BTreeMap<String, f64>,
    tol: f64,
    label: &str,
) {
    let only_old: Vec<&String> = old.keys().filter(|k| !new.contains_key(*k)).collect();
    let only_new: Vec<&String> = new.keys().filter(|k| !old.contains_key(*k)).collect();
    assert!(
        only_old.is_empty() && only_new.is_empty(),
        "[{label}] support key sets differ: {} only in old (e.g. {:?}), {} only in new (e.g. {:?})",
        only_old.len(),
        only_old.first().map(|k| (k, old[*k])),
        only_new.len(),
        only_new.first().map(|k| (k, new[*k])),
    );
    for (k, o) in old {
        let n = new[k];
        assert!(
            (o - n).abs() <= tol,
            "[{label}] coefficient at {k} differs: old {o} vs new {n} (tol {tol})"
        );
    }
    assert!(!old.is_empty(), "[{label}] final support is empty");
}

// ---------------------------------------------------------------------------
// Workload 1 — trotter-tfim-noisy: golden master + truncation fidelity.
// ---------------------------------------------------------------------------

/// Old's `tests/trotter.rs` acceptance bar, reproduced on BOTH engines.
///
/// * (1) the exact (`NoStrategy`/`NoPolicy`) `N = 4`, `STEPS = 10` expectation
///   equals the frozen `GOLDEN` to `1e-9`;
/// * (2) the `CoefficientThreshold(1e-6)` run stays within `1e-6` of exact.
#[test]
fn trotter_golden_master_and_truncation_fidelity_hold_on_both_engines() {
    const N: usize = 4;
    const STEPS: usize = 10;
    const THETA_X: f64 = 0.1; // dt·h, dt = 0.1, h = 1
    const THETA_ZZ: f64 = 0.0125; // dt·J, J = 1/8
    const NOISE: [f64; 3] = [2.5e-5; 3]; // 1e-4 / 4 per channel
    /// Old's frozen constant (`ppvm-pauli-sum/tests/trotter.rs:104`).
    const GOLDEN: f64 = 2.1610566562692544;

    // --- Exact runs (no truncation at all). ---------------------------------
    let mut old_exact: OldExactSum = OldPauliSum::builder().n_qubits(N).build();
    seed_old!(old_exact, N);
    trotter_evolve!(old_exact, N, STEPS, THETA_X, THETA_ZZ, NOISE);
    let old_exact_val = old_zero_expectation!(old_exact);

    let mut new_exact = NewExactSum::with_policy(N, NoPolicy);
    seed_new!(new_exact, N);
    trotter_evolve!(new_exact, N, STEPS, THETA_X, THETA_ZZ, NOISE);
    let new_exact_val = new_zero_expectation!(new_exact);

    // (1) Golden master. Asserted on old too: it pins that the `cnot; rz; cnot`
    // decomposition used on both sides reproduces old's native-`rzz` constant,
    // so the new engine is being held to old's actual number.
    assert!(
        (old_exact_val - GOLDEN).abs() < 1e-9,
        "OLD exact Trotter expectation {old_exact_val} drifted from golden {GOLDEN}"
    );
    assert!(
        (new_exact_val - GOLDEN).abs() < 1e-9,
        "NEW exact Trotter expectation {new_exact_val} drifted from golden {GOLDEN}"
    );

    // The two exact runs must also agree with each other far below the golden
    // tolerance (same circuit, same arithmetic, only merge order differs).
    assert_supports_identical(
        &old_support!(old_exact),
        &new_support!(new_exact),
        1e-12,
        "trotter N=4 exact",
    );

    // --- Truncated runs (coefficient floor only, `1e-6`). -------------------
    let mut old_approx: OldThreshSum = OldPauliSum::builder()
        .n_qubits(N)
        .strategy(OldCoeffThreshold(THRESHOLD))
        .build();
    seed_old!(old_approx, N);
    trotter_evolve!(old_approx, N, STEPS, THETA_X, THETA_ZZ, NOISE);
    let old_approx_val = old_zero_expectation!(old_approx);

    let mut new_approx = NewThreshSum::with_policy(
        N,
        NewCoeffThreshold {
            threshold: THRESHOLD,
        },
    );
    seed_new!(new_approx, N);
    trotter_evolve!(new_approx, N, STEPS, THETA_X, THETA_ZZ, NOISE);
    let new_approx_val = new_zero_expectation!(new_approx);

    // (2) Truncation fidelity, on both engines. An engine that pruned at
    // insertion (rather than deferring truncation to the caller) drifts ~2.2e-6.
    let old_drift = (old_exact_val - old_approx_val).abs();
    let new_drift = (new_exact_val - new_approx_val).abs();
    assert!(
        old_drift < 1e-6,
        "OLD truncated result drifted from exact by {old_drift:e}"
    );
    assert!(
        new_drift < 1e-6,
        "NEW truncated result drifted from exact by {new_drift:e} — an engine that \
         prunes at insertion drifts ~2.2e-6"
    );
    // The two engines must truncate identically, so their drifts coincide to
    // round-off — a much sharper statement than each being under `1e-6`.
    assert!(
        (old_drift - new_drift).abs() < 1e-12,
        "truncation drift differs: old {old_drift:e} vs new {new_drift:e}"
    );
    assert_supports_identical(
        &old_support!(old_approx),
        &new_support!(new_approx),
        1e-12,
        "trotter N=4 truncated",
    );
}

/// The headline workload at its benchmarked width: `n = 12`, 10 Trotter steps,
/// `CombinedStrategy/Policy(CoefficientThreshold(1e-6), MaxPauliWeight(MAX))`,
/// capacity `n²`. Strict whole-support diff (exact key-set equality) plus the
/// zero-state expectation.
#[test]
fn trotter_n12_support_matches_old_exactly() {
    let n = 12usize;
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let steps = ((1.0 / h) / dt) as usize;
    let theta_x = dt * h;
    let theta_zz = dt * (1.0 / 8.0 * h);
    let noise = [1e-4 / 4.0; 3];

    let mut old: OldCombinedSum = OldPauliSum::builder()
        .n_qubits(n)
        .strategy(CombinedStrategy(
            OldCoeffThreshold(THRESHOLD),
            OldMaxWeight(usize::MAX),
        ))
        .capacity(n.pow(2))
        .build();
    seed_old!(old, n);

    let mut new = NewCombinedSum::with_capacity(
        n,
        CombinedPolicy(
            NewCoeffThreshold {
                threshold: THRESHOLD,
            },
            NewMaxWeight(usize::MAX),
        ),
        n.pow(2),
    );
    seed_new!(new, n);

    trotter_evolve!(old, n, steps, theta_x, theta_zz, noise);
    trotter_evolve!(new, n, steps, theta_x, theta_zz, noise);

    let os = old_support!(old);
    let ns = new_support!(new);
    // Well below the `1e-6` floor, so no term can be reclassified kept/dropped.
    assert_supports_identical(&os, &ns, 1e-9, "trotter n=12");
    assert!(
        os.len() > 100,
        "the n=12 workload should grow a real support, got {}",
        os.len()
    );

    let ov = old_zero_expectation!(old);
    let nv = new_zero_expectation!(new);
    assert!(
        (ov - nv).abs() <= 1e-9 * ov.abs().max(1.0),
        "zero-state expectation differs: old {ov} vs new {nv}"
    );
}

// ---------------------------------------------------------------------------
// Workload 2 — trotter-qubit-scaling sweep (coefficient floor alone).
// ---------------------------------------------------------------------------

/// The same circuit swept over qubit count under `CoefficientThreshold(1e-6)`
/// **alone** (no `CombinedStrategy`), capacity `n²`, driven at `J = 1.0` so the
/// support runs large. The bar is *exact* support-size equality at every `n` —
/// with a coefficient floor the surviving count is a sensitive function of every
/// accumulated coefficient, so an engine that merges or truncates even slightly
/// differently lands on a different integer.
///
/// The sweep is `n = 2, 6, 10, 14, 18` (old's `(2..65).step_by(4)` shape, capped
/// so a debug-profile test stays quick; the bench sweeps the wide end).
#[test]
fn trotter_qubit_sweep_support_and_expectation_match_old() {
    let h = 1.0_f64;
    let dt = 0.1 / h;
    let steps = 5usize;
    let theta_x = dt * h;
    let theta_zz = dt * 1.0; // J = 1.0: drive the support large
    let noise = [1e-4 / 4.0; 3];

    for n in (2..19usize).step_by(4) {
        let mut old: OldThreshSum = OldPauliSum::builder()
            .n_qubits(n)
            .strategy(OldCoeffThreshold(THRESHOLD))
            .capacity(n.pow(2))
            .build();
        seed_old!(old, n);

        let mut new = NewThreshSum::with_capacity(
            n,
            NewCoeffThreshold {
                threshold: THRESHOLD,
            },
            n.pow(2),
        );
        seed_new!(new, n);

        trotter_evolve!(old, n, steps, theta_x, theta_zz, noise);
        trotter_evolve!(new, n, steps, theta_x, theta_zz, noise);

        assert_eq!(
            old.data().len(),
            new.len(),
            "[sweep n={n}] final |support| differs: old {} vs new {}",
            old.data().len(),
            new.len()
        );
        assert_supports_identical(
            &old_support!(old),
            &new_support!(new),
            1e-9,
            &format!("sweep n={n}"),
        );

        let ov = old_zero_expectation!(old);
        let nv = new_zero_expectation!(new);
        assert!(
            (ov - nv).abs() <= 1e-9 * ov.abs().max(1.0),
            "[sweep n={n}] expectation differs: old {ov} vs new {nv}"
        );
    }
}

// ---------------------------------------------------------------------------
// Workload 3 — deep random circuit with NO truncation (pure fan-out growth).
// ---------------------------------------------------------------------------

/// One `depth`-layer rotation+entangler circuit: `rz(1.1); ry(2.1); rz(1.1)` on
/// every qubit, then a `cnot` ring, repeated `depth` times and closed with a
/// final rotation layer (old's `benches/random-circuit.rs` shape). No truncation
/// anywhere — the support grows monotonically by fan-out.
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
fn untruncated_deep_circuit_matches_old_state_and_norm() {
    // `n = 8`, `depth = 2`: fan-out with no truncation to shrink it, so this is a
    // strict whole-state diff over a support in the thousands. (The old bench's
    // committed params are `n = 4, depth = 2`; its gate-bench prologue uses
    // `n = 12, depth = 2`. Debug-profile runtime caps this at 8.)
    let n = 8usize;
    let depth = 2usize;

    let zz: String = (0..n)
        .map(|i| if i < 2 { 'Z' } else { 'I' })
        .collect::<String>();

    let mut old: OldExactSum = OldPauliSum::builder().n_qubits(n).build();
    old += (zz.as_str(), 1.0);
    let mut new = NewExactSum::with_policy(n, NoPolicy);
    new += (NewKey::from(zz.as_str()), 1.0);

    random_circuit_evolve!(old, n, depth);
    random_circuit_evolve!(new, n, depth);

    let os = old_support!(old);
    let ns = new_support!(new);
    assert_supports_identical(&os, &ns, 1e-9, "untruncated deep circuit");
    assert!(
        os.len() > 500,
        "untruncated fan-out should reach a large support, got {}",
        os.len()
    );

    // The ℓ² norm — `overlap(f, f)` — is invariant under key ordering, so it
    // catches a dropped or duplicated term that a per-key loop could miss.
    let old_norm: f64 = os.values().map(|c| c * c).sum();
    let new_norm = new.overlap(&new);
    assert!(
        (old_norm - new_norm).abs() <= 1e-9 * old_norm.abs().max(1.0),
        "ℓ² norm differs: old {old_norm} vs new {new_norm}"
    );
}

/// The seeded random *gate stream* (the conformance generator) replayed with **no
/// truncation**: a heterogeneous `H`/`S`/`CNOT`/`Rx`/`Rz` order rather than the
/// structured layer/ring shape above. Kept short and narrow because untruncated
/// fan-out is exponential in the rotation count.
#[test]
fn untruncated_random_gate_stream_matches_old() {
    use ppvm_conformance_2::{GateOp, random_circuit};

    for &seed in &[1u64, 42, 777] {
        let mut rng: StdRng = seeded_rng(seed);
        let n = 6usize;
        let circuit = random_circuit(&mut rng, n, 40);

        let zz: String = (0..n).map(|i| if i < 2 { 'Z' } else { 'I' }).collect();
        let mut old: OldExactSum = OldPauliSum::builder().n_qubits(n).build();
        old += (zz.as_str(), 1.0);
        let mut new = NewExactSum::with_policy(n, NoPolicy);
        new += (NewKey::from(zz.as_str()), 1.0);

        for &op in &circuit {
            match op {
                GateOp::H(q) => {
                    old.h(q);
                    new.h(q);
                }
                GateOp::S(q) => {
                    old.s(q);
                    new.s(q);
                }
                GateOp::Cnot(c, t) => {
                    old.cnot(c, t);
                    new.cnot(c, t);
                }
                GateOp::Rx(q, th) => {
                    old.rx(q, th);
                    new.rx(q, th);
                }
                GateOp::Rz(q, th) => {
                    old.rz(q, th);
                    new.rz(q, th);
                }
            }
        }

        let os = old_support!(old);
        let ns = new_support!(new);
        assert_supports_identical(&os, &ns, 1e-9, &format!("untruncated stream seed={seed}"));

        let old_norm: f64 = os.values().map(|c| c * c).sum();
        let new_norm = new.overlap(&new);
        assert!(
            (old_norm - new_norm).abs() <= 1e-9 * old_norm.abs().max(1.0),
            "[seed={seed}] ℓ² norm differs: old {old_norm} vs new {new_norm}"
        );
    }
}
