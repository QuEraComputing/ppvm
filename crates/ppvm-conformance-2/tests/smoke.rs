// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phase 0 smoke tests: prove the harness mechanics work end to end — the
//! generators are deterministic under a seed, and they feed the old crate. The
//! real old-vs-new differential suites arrive with the `-2` data types.

use ppvm_conformance_2::{
    GateOp, assert_close, old_word_from_str, random_circuit, random_pauli_string, seeded_rng,
};

#[test]
fn generators_are_deterministic() {
    // Same seed => identical output, so old and new backends can be driven by
    // the same randomness in a differential test.
    let mut a = seeded_rng(42);
    let mut b = seeded_rng(42);
    assert_eq!(
        random_pauli_string(&mut a, 16),
        random_pauli_string(&mut b, 16)
    );

    let mut c = seeded_rng(7);
    let mut d = seeded_rng(7);
    let cc: Vec<GateOp> = random_circuit(&mut c, 5, 128);
    let dd: Vec<GateOp> = random_circuit(&mut d, 5, 128);
    assert_eq!(cc, dd);
    assert_eq!(cc.len(), 128);
}

#[test]
fn different_seeds_differ() {
    let s1 = random_pauli_string(&mut seeded_rng(1), 16);
    let s2 = random_pauli_string(&mut seeded_rng(2), 16);
    assert_ne!(s1, s2);
}

#[test]
fn feeds_the_old_crate() {
    // The generated string builds a real old `PauliWord`; construction is
    // deterministic. Phase 2 adds the `-2` twin and diffs against it.
    let s = random_pauli_string(&mut seeded_rng(1), 8);
    let w1 = old_word_from_str(&s);
    let w2 = old_word_from_str(&s);
    assert_eq!(w1, w2);
}

#[test]
fn assert_close_within_tolerance() {
    assert_close(1.0, 1.0 + 1e-12, 1e-9);
}
