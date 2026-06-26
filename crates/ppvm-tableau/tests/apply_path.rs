// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Differential lock for the two-qubit-rotation (`rotate_2`) coefficient path.
//!
//! `rotate_2` (RXX/RYY/RZZ) is the *only* caller of
//! `compute_coefficients_after_pauli_apply`, i.e. the "apply" coefficient
//! accumulation in `data.rs`. That accumulation relabels every branch index by
//! a fixed `idx ^ stab_anticomm_bits`, which is a bijection — so the keys never
//! collide and the coalescing container can never actually merge two entries.
//!
//! This test pins the measured-bit record of a branchy RXX/RYY/RZZ brickwork
//! over many seeds to an FNV-1a digest, so that swapping the apply
//! accumulation's storage strategy (hash coalesce → direct relabel) is proven
//! to leave every measurement outcome bit-identical.

use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ByteFxHashF64<8>, usize>;

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fnv1a_update(mut h: u64, byte: u8) -> u64 {
    h ^= byte as u64;
    h = h.wrapping_mul(FNV_PRIME);
    h
}

/// Deterministic branchy two-qubit-rotation circuit on `n` qubits.
///
/// Brickwork layers of RXX / RYY / RZZ at non-Clifford angles, interleaved with
/// Hadamards, so the coefficient vector genuinely branches and every
/// `rotate_2` exercises the apply path on a non-trivial superposition.
fn build_rot2_brickwork(n: usize, layers: usize) -> Tab {
    let mut tab: Tab = GeneralizedTableau::new_with_seed(n, 1e-10, 1);
    for q in (0..n).step_by(2) {
        tab.h(q);
    }
    for layer in 0..layers {
        // even brickwork pairs
        for a in (0..n.saturating_sub(1)).step_by(2) {
            tab.rxx(a, a + 1, 0.3 * std::f64::consts::PI);
            tab.ryy(a, a + 1, 0.4 * std::f64::consts::PI);
        }
        // odd brickwork pairs
        for a in (1..n.saturating_sub(1)).step_by(2) {
            tab.rzz(a, a + 1, 0.25 * std::f64::consts::PI);
            tab.rxx(a, a + 1, 0.15 * std::f64::consts::PI);
        }
        if layer % 2 == 0 {
            for q in (1..n).step_by(2) {
                tab.h(q);
            }
        }
    }
    tab
}

/// Fork `tab` over `seeds` independent RNG streams, measure every qubit, and
/// fold the full outcome record into an FNV-1a digest.
fn measure_record_digest(tab: &Tab, n: usize, seeds: u64) -> u64 {
    let mut h = FNV_OFFSET;
    for seed in 0..seeds {
        let mut forked = tab.fork(Some(seed));
        for q in 0..n {
            let bit = forked.measure(q).expect("no lost qubits in this circuit");
            h = fnv1a_update(h, bit as u8);
        }
    }
    h
}

#[test]
fn rot2_apply_path_measurement_digest_is_stable() {
    let n = 8;
    let tab = build_rot2_brickwork(n, 3);
    // The circuit must actually branch, or it wouldn't exercise the apply path.
    assert!(
        tab.coefficients.len() > 8,
        "expected a branchy superposition, got {} coefficients",
        tab.coefficients.len()
    );

    let digest = measure_record_digest(&tab, n, 256);
    println!("rot2_apply_path digest = {digest:#018x}");

    // Golden digest captured on the hash-coalesce apply path (pre-refactor).
    // The direct-relabel apply path must reproduce it bit-for-bit.
    assert_eq!(
        digest, 0x2401_e08e_70e6_ecc8,
        "measurement record changed — apply-path refactor is not behaviour-preserving"
    );
}
