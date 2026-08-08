// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Test-only conformance & benchmark harness for the `ppvm-*-2` refactor.
//!
//! It provides the pieces every differential test and comparative benchmark
//! reuses: a seeded RNG (so the old and new backends observe *identical*
//! randomness), generators for random Pauli words and Clifford+rotation
//! circuits emitted as replayable data, and coefficient-comparison helpers.
//!
//! At Phase 0 (scaffolding) only the generators and a single old-crate
//! constructor exist; the per-crate differential suites are added as the `-2`
//! crates come online (a `PauliWord` twin diff in Phase 2, `Sum` in Phase 3, …).

use ppvm_pauli_word::prelude::PauliWord;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// A deterministic RNG. Seeding old and new backends from the same value makes
/// their randomness identical, which is what lets a differential test compare
/// them term for term.
pub fn seeded_rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

/// The four single-qubit Pauli letters, in `Word::Site` order.
pub const PAULIS: [char; 4] = ['I', 'X', 'Y', 'Z'];

/// A random Pauli word of `n` qubits as a string like `"XIZY"`. Reusable to
/// build either an old [`PauliWord`] (via [`old_word_from_str`]) or a future
/// `-2` word from the *same* string, so the two can be diffed.
pub fn random_pauli_string(rng: &mut StdRng, n: usize) -> String {
    (0..n)
        .map(|_| PAULIS[rng.random_range(0..4usize)])
        .collect()
}

/// The five lossy Pauli letters (`I`/`X`/`Y`/`Z` and the loss symbol `L`), in the
/// symbol order both the old `LossyPauliWord` and the new `-2` one parse.
pub const LOSSY_LETTERS: [char; 5] = ['I', 'X', 'Y', 'Z', 'L'];

/// A random **lossy** Pauli word of `n` qubits as a string like `"XILZL"`: each
/// site is independently one of `I`/`X`/`Y`/`Z`/`L` (so ≈1/5 of sites are `Lost`),
/// giving a spread that includes fully-present, mixed, and fully-lost words.
///
/// Reusable to build the old `ppvm-pauli-word::LossyPauliWord` and the new
/// `ppvm-lossy-pauli-word-2::LossyPauliWord` from the *same* string so the two can
/// be diffed site-for-site, including their loss planes.
pub fn random_lossy_pauli_string(rng: &mut StdRng, n: usize) -> String {
    (0..n)
        .map(|_| LOSSY_LETTERS[rng.random_range(0..5usize)])
        .collect()
}

/// A replayable gate operation: emitted once as data, then applied to any
/// backend (old or new) so both see the same circuit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GateOp {
    /// Hadamard on the given qubit.
    H(usize),
    /// Phase gate `S` on the given qubit.
    S(usize),
    /// CNOT with `(control, target)`.
    Cnot(usize, usize),
    /// X-rotation by the given angle on the given qubit.
    Rx(usize, f64),
    /// Z-rotation by the given angle on the given qubit.
    Rz(usize, f64),
}

/// A random Clifford+rotation circuit on `n_qubits`, `len` gates long. Returned
/// as plain data so the identical circuit can drive an old and a new backend.
pub fn random_circuit(rng: &mut StdRng, n_qubits: usize, len: usize) -> Vec<GateOp> {
    assert!(n_qubits >= 1, "a circuit needs at least one qubit");
    let angle = |rng: &mut StdRng| rng.random_range(-std::f64::consts::PI..std::f64::consts::PI);
    (0..len)
        .map(|_| {
            let q = rng.random_range(0..n_qubits);
            match rng.random_range(0..5usize) {
                0 => GateOp::H(q),
                1 => GateOp::S(q),
                2 if n_qubits > 1 => {
                    let mut t = rng.random_range(0..n_qubits);
                    while t == q {
                        t = rng.random_range(0..n_qubits);
                    }
                    GateOp::Cnot(q, t)
                }
                2 => GateOp::H(q), // single-qubit register: no CNOT target
                3 => GateOp::Rx(q, angle(rng)),
                _ => GateOp::Rz(q, angle(rng)),
            }
        })
        .collect()
}

/// The old-crate word type used by the harness: byte-backed, 16-qubit capacity,
/// default (`fxhash`) hasher — matching the shapes the current crates test with.
pub type OldWord = PauliWord<[u8; 2]>;

/// Build an old-crate [`PauliWord`] from a generated Pauli string. The `-2` twin
/// is diffed against this starting in Phase 2.
pub fn old_word_from_str(s: &str) -> OldWord {
    OldWord::from(s)
}

/// Assert two `f64` coefficients agree within tolerance. Differential tests on
/// floating-point backends compare within `tol` rather than bit-for-bit.
#[track_caller]
pub fn assert_close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "coefficients differ: {a} vs {b} (tol {tol})"
    );
}

// ===========================================================================
// Phase-3 `Sum` differential harness — build matched OLD and NEW `PauliSum<f64>`
// from the *same* `(pauli_string, coeff)` term list and compare observables.
// ===========================================================================

use ppvm_pauli_sum::config::fxhash::ByteF64;
use ppvm_pauli_sum::strategy::CoefficientThreshold;
use ppvm_pauli_sum::sum::PauliSum as OldPauliSum;
use ppvm_pauli_sum_2::{PauliSum as NewPauliSum, PauliWord as NewPauliWord};

/// The OLD reference `PauliSum<f64>`: `[u8; 8]` (64-qubit-capacity) storage,
/// `FxHash`, and a [`CoefficientThreshold`] strategy so its `truncate()` realizes
/// the NEW crate's structural `reduce` (drop zero-coefficient keys). The default
/// threshold is `1e-12`; exact cancellations (`1.0 + (-1.0) = 0.0`) fall under it,
/// so the two crates drop the same keys.
pub type OldSum = OldPauliSum<ByteF64<8, CoefficientThreshold>>;

/// The NEW `PauliSum<f64, NoPolicy>` under test — the **shipped default**
/// (`u64`-backed [`PauliWord`], `IdentityBuildHasher`). The differential and
/// Lean-oracle suites run against this so they exercise the configuration that
/// actually ships. (The perf benchmark separately pins its *own* `[u8; 8]`-backed
/// variant to storage-match [`OldSum`] for a fair engine-to-engine ratio — that is
/// a bench-local concern, since correctness is storage-independent.)
pub type NewSum = NewPauliSum;

/// The NEW crate's key type (`PauliWord<u64>`, 64-qubit-capacity), re-exported so
/// tests can name keys for `contains`/`get`.
pub type NewKey = NewPauliWord;

/// Build the OLD reference sum on `n_qubits` from `(string, coeff)` terms.
///
/// Colliding keys are combined by the old `+=` (`add_assign`); zero-coefficient
/// keys are **not** dropped here (call [`reduce_old`] to realize `reduce`).
pub fn build_old_sum(n_qubits: usize, terms: &[(String, f64)]) -> OldSum {
    let mut s: OldSum = OldPauliSum::builder().n_qubits(n_qubits).build();
    for (w, c) in terms {
        s += (w.as_str(), *c);
    }
    s
}

/// Build the NEW sum on `n_qubits` from the *same* `(string, coeff)` terms.
///
/// `from_terms` runs `accumulate_batch` only: colliding keys are combined and a
/// zero (or exactly-cancelling) coefficient is **kept**, matching old's `+=`.
/// Canonicalization is caller-driven on both sides ([`reduce_old`] /
/// `Sum::reduce`).
pub fn build_new_sum(n_qubits: usize, terms: &[(String, f64)]) -> NewSum {
    NewSum::from_terms(
        n_qubits,
        terms
            .iter()
            .map(|(w, c)| (NewPauliWord::from(w.as_str()), *c)),
    )
}

/// Realize the NEW crate's structural `reduce` on the OLD sum by running its
/// coefficient-threshold `truncate()` (drops the exact-zero keys a cancellation
/// leaves behind).
pub fn reduce_old(s: &mut OldSum) {
    s.truncate();
}

/// The OLD sum's support as a sorted `(canonical_pauli_string, coeff)` vector.
pub fn old_support(s: &OldSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.data().iter().map(|(k, c)| (k.to_string(), *c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// The NEW sum's support as a sorted `(canonical_pauli_string, coeff)` vector.
pub fn new_support(s: &NewSum) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = s.iter().map(|(k, c)| (k.to_string(), c)).collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// Assert the OLD and NEW supports agree as sorted `(string, coeff)` sets within
/// `tol`: same keys in the same canonical order, coefficients close.
#[track_caller]
pub fn assert_supports_match(old: &OldSum, new: &NewSum, tol: f64) {
    let os = old_support(old);
    let ns = new_support(new);
    assert_eq!(
        os.len(),
        ns.len(),
        "support size differs: old {} vs new {}\nold={os:?}\nnew={ns:?}",
        os.len(),
        ns.len()
    );
    for (o, n) in os.iter().zip(ns.iter()) {
        assert_eq!(o.0, n.0, "support key differs: old {} vs new {}", o.0, n.0);
        assert_close(o.1, n.1, tol);
    }
}

/// A random distinct-key term list of `count` `(pauli_string, coeff)` pairs on
/// `n` qubits, coefficients in `[-2, 2)`. Keys may repeat (both crates merge them
/// identically); coefficients avoid the near-zero band so no term accidentally
/// falls under the OLD threshold.
pub fn random_terms(rng: &mut StdRng, n: usize, count: usize) -> Vec<(String, f64)> {
    (0..count)
        .map(|_| {
            let w = random_pauli_string(rng, n);
            // Coefficient in [-2, -0.25] ∪ [0.25, 2): O(1), never near the 1e-12
            // reduce threshold.
            let mag = rng.random_range(0.25..2.0);
            let c = if rng.random_range(0..2usize) == 0 {
                mag
            } else {
                -mag
            };
            (w, c)
        })
        .collect()
}

pub mod mixture;
/// Phase-4 tableau differential harness: matched OLD/NEW engines behind one
/// [`Driver`](tableau::Driver) trait, plus the integration-baseline workloads.
pub mod tableau;

/// Phase-5 symbolic-coefficient differential harness: matched OLD/NEW
/// `Term`-coefficient sums and the `sym.*` integration-baseline workloads.
pub mod sym;
