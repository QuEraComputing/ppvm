// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Pauli-word width is a const-generic parameter, so a `LindbladSpec` can be
//! instantiated wider than the historical 128-qubit ceiling.
//!
//! The central check is a *padding invariance*: the same physical problem,
//! embedded in words of different widths, must produce bit-identical numbers.
//! A width bug (mask truncation, a stray `W_CHUNKS`, a wrong support index)
//! breaks that immediately.

use num::Complex;
use ppvm_lindblad::{JumpInput, LindbladSpec, PcStepConfig, max_qubits, word_from_codes};

/// `H = J Σ_{i<i+1} Z_iZ_{i+1} + h Σ_i X_i` on the first `n_active` qubits of
/// an `n_total`-qubit register, plus per-site dephasing. The remaining qubits
/// are spectators: identity everywhere.
fn model(n_total: usize, n_active: usize) -> (Vec<(String, f64)>, Vec<JumpInput>) {
    let pad = |sites: &[(usize, char)]| -> String {
        let mut s = vec!['I'; n_total];
        for &(q, c) in sites {
            s[q] = c;
        }
        s.into_iter().collect()
    };
    let mut h = Vec::new();
    for i in 0..n_active - 1 {
        h.push((pad(&[(i, 'Z'), (i + 1, 'Z')]), 0.7));
    }
    for i in 0..n_active {
        h.push((pad(&[(i, 'X')]), 1.3));
    }
    let jumps = (0..n_active)
        .map(|i| JumpInput {
            lincomb: vec![(pad(&[(i, 'Z')]), Complex::new(1.0, 0.0))],
            rate: 0.05,
        })
        .collect();
    (h, jumps)
}

/// Evolve `Z_0 Z_1` for a few steps and return the coefficient sum over
/// `{I, Z}`-only strings — i.e. the expectation on the all-`Z = -1` product
/// state, up to the sign convention (identical across widths, which is all
/// this test needs).
fn evolve<const C: usize>(n_total: usize, n_active: usize, steps: usize) -> f64 {
    let (h, jumps) = model(n_total, n_active);
    let spec = LindbladSpec::<C>::new(n_total, &h, &jumps).unwrap();

    let mut codes = vec![0u8; n_total];
    codes[0] = 2; // Z
    codes[1] = 2; // Z
    let mut basis = vec![word_from_codes::<C>(&codes).unwrap()];
    let mut coeffs = vec![1.0f64];

    let cfg = PcStepConfig {
        max_basis: 20_000,
        admit_basis: Some(60_000),
        drop_tol: 0.0,
        tau_add: None,
        num_threads: Some(1),
    };
    for _ in 0..steps {
        spec.pc_step(&mut basis, &mut coeffs, 0.05, &[], &cfg)
            .unwrap();
    }

    let mut out = vec![0u8; n_total];
    let mut acc = 0.0;
    for (w, c) in basis.iter().zip(&coeffs) {
        ppvm_lindblad::codes_from_word(w, &mut out);
        if out.iter().all(|&b| b == 0 || b == 2) {
            let nz = out.iter().filter(|&&b| b == 2).count();
            acc += if nz % 2 == 0 { *c } else { -*c };
        }
    }
    acc
}

#[test]
fn capacity_scales_with_chunk_count() {
    assert_eq!(max_qubits::<2>(), 128);
    assert_eq!(max_qubits::<4>(), 256);
    assert_eq!(max_qubits::<8>(), 512);
}

#[test]
fn wide_words_exceed_the_old_128_qubit_ceiling() {
    // The exact case that used to fail with "supports n_qubits <= 128".
    let (h, jumps) = model(256, 4);
    let spec = LindbladSpec::<4>::new(256, &h, &jumps).unwrap();
    assert_eq!(spec.n_qubits(), 256);

    let (h, jumps) = model(512, 4);
    let spec = LindbladSpec::<8>::new(512, &h, &jumps).unwrap();
    assert_eq!(spec.n_qubits(), 512);
}

#[test]
fn too_many_qubits_for_the_width_is_rejected() {
    let (h, jumps) = model(8, 4);
    // 8 qubits of content declared as a 200-qubit register: C = 2 holds 128.
    let padded: Vec<(String, f64)> = h
        .iter()
        .map(|(s, c)| (format!("{s}{}", "I".repeat(192)), *c))
        .collect();
    let jumps: Vec<JumpInput> = jumps
        .iter()
        .map(|j| JumpInput {
            lincomb: j
                .lincomb
                .iter()
                .map(|(s, c)| (format!("{s}{}", "I".repeat(192)), *c))
                .collect(),
            rate: j.rate,
        })
        .collect();
    assert!(LindbladSpec::<2>::new(200, &padded, &jumps).is_err());
    assert!(LindbladSpec::<4>::new(200, &padded, &jumps).is_ok());
}

#[test]
fn padding_into_a_wider_word_changes_nothing() {
    // Identical 6-qubit physics, embedded in 64-, 200- and 400-qubit
    // registers backed by 2-, 4- and 8-chunk words.
    let narrow = evolve::<2>(64, 6, 6);
    let wide = evolve::<4>(200, 6, 6);
    let widest = evolve::<8>(400, 6, 6);

    assert!(narrow.abs() > 1e-6, "test observable is trivially zero");
    assert_eq!(
        narrow.to_bits(),
        wide.to_bits(),
        "C=2 gave {narrow}, C=4 gave {wide}"
    );
    assert_eq!(
        narrow.to_bits(),
        widest.to_bits(),
        "C=2 gave {narrow}, C=8 gave {widest}"
    );
}

#[test]
fn same_width_different_register_size_agrees() {
    // Within one width, the spectator qubits must not touch the answer.
    assert_eq!(
        evolve::<4>(130, 6, 5).to_bits(),
        evolve::<4>(256, 6, 5).to_bits()
    );
}
