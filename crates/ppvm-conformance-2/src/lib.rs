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
