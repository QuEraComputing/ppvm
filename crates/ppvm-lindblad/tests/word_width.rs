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
use ppvm_lindblad::{
    CHUNK_BITS, JumpInput, LindbladSpec, MAX_SUPPORTED_QUBITS, PcStepConfig, Sector, W_CHUNKS,
    WIDTHS, Word, chunks_for, codes_from_word, max_qubits, word_from_codes,
};
use ppvm_pauli_sum::symmetry::{TranslationGroup, canonicalize_pauli_sum_complex};

const W128: usize = WIDTHS[0];
const W256: usize = WIDTHS[1];
const W512: usize = WIDTHS[2];

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

/// A decay-type Kossakowski dissipator `A_0 = σ⁻_0`, `A_1 = σ⁻_1` with a
/// non-diagonal pair matrix, on the first two qubits of `n_total`.
fn add_kossakowski<const C: usize>(spec: &mut LindbladSpec<C>, n_total: usize) {
    let op = |q: usize| -> Vec<(String, Complex<f64>)> {
        let mut x = vec!['I'; n_total];
        let mut y = vec!['I'; n_total];
        x[q] = 'X';
        y[q] = 'Y';
        vec![
            (x.into_iter().collect(), Complex::new(0.5, 0.0)),
            (y.into_iter().collect(), Complex::new(0.0, -0.5)),
        ]
    };
    let k = vec![
        vec![Complex::new(0.3, 0.0), Complex::new(0.1, 0.05)],
        vec![Complex::new(0.1, -0.05), Complex::new(0.2, 0.0)],
    ];
    spec.add_kossakowski(&[op(0), op(1)], &k).unwrap();
}

/// Evolve `Z_0 Z_1` for a few steps and return the coefficient sum over
/// `{I, Z}`-only strings — i.e. the expectation on the all-`Z = -1` product
/// state, up to the sign convention (identical across widths, which is all
/// this test needs).
fn evolve<const C: usize>(n_total: usize, n_active: usize, steps: usize, koss: bool) -> f64 {
    let (h, jumps) = model(n_total, n_active);
    let mut spec = LindbladSpec::<C>::new(n_total, &h, &jumps).unwrap();
    if koss {
        add_kossakowski(&mut spec, n_total);
    }

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
        codes_from_word(w, &mut out);
        if out.iter().all(|&b| b == 0 || b == 2) {
            let nz = out.iter().filter(|&&b| b == 2).count();
            acc += if nz % 2 == 0 { *c } else { -*c };
        }
    }
    acc
}

/// Momentum-sector evolution of `Σ_q X_q` under a translation-invariant ring
/// `H = Σ (X X + Y Y + ½ Z Z) + 0.3 Σ Z` of `n` sites, untruncated. Returns
/// the rep coefficients sorted by word.
fn evolve_orbit<const C: usize>(n: usize, k: i32, steps: usize) -> Vec<(Vec<u8>, Complex<f64>)> {
    let word = |ops: &[(usize, char)]| -> String {
        let mut s = vec!['I'; n];
        for &(q, c) in ops {
            s[q] = c;
        }
        s.into_iter().collect()
    };
    let mut h = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        h.push((word(&[(i, 'X'), (j, 'X')]), 1.0));
        h.push((word(&[(i, 'Y'), (j, 'Y')]), 1.0));
        h.push((word(&[(i, 'Z'), (j, 'Z')]), 0.5));
        h.push((word(&[(i, 'Z')]), 0.3));
    }
    let spec = LindbladSpec::<C>::new(n, &h, &[]).unwrap();
    let group = TranslationGroup::chain_1d(n);
    let k_modes = [k];

    let mut basis: Vec<Word<C>> = (0..n)
        .map(|q| {
            let mut codes = vec![0u8; n];
            codes[q] = 1; // X
            word_from_codes::<C>(&codes).unwrap()
        })
        .collect();
    let mut coeffs: Vec<Complex<f64>> = (0..n)
        .map(|q| {
            Complex::from_polar(
                1.0,
                -2.0 * std::f64::consts::PI * (k as f64) * q as f64 / n as f64,
            )
        })
        .collect();
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &k_modes);

    let sector = Sector::new(&group, &k_modes);
    let cfg = PcStepConfig {
        max_basis: usize::MAX,
        admit_basis: None,
        drop_tol: 0.0,
        tau_add: None,
        num_threads: Some(1),
    };
    for _ in 0..steps {
        spec.pc_step_orbit_rep(&mut basis, &mut coeffs, 0.1, &[], &sector, &cfg)
            .unwrap();
    }
    let mut out: Vec<(Vec<u8>, Complex<f64>)> = basis
        .iter()
        .zip(&coeffs)
        .map(|(w, c)| {
            let mut codes = vec![0u8; n];
            codes_from_word(w, &mut codes);
            (codes, *c)
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[test]
fn capacity_scales_with_chunk_count() {
    assert_eq!(max_qubits::<W128>(), 128);
    assert_eq!(max_qubits::<W256>(), 256);
    assert_eq!(max_qubits::<W512>(), 512);
    assert_eq!(W_CHUNKS, W128);
    assert_eq!(W128 * CHUNK_BITS, 128);
    assert_eq!(chunks_for(1), Some(W128));
    assert_eq!(chunks_for(128), Some(W128));
    assert_eq!(chunks_for(129), Some(W256));
    assert_eq!(chunks_for(130), Some(W256));
    assert_eq!(chunks_for(512), Some(W512));
    assert_eq!(chunks_for(MAX_SUPPORTED_QUBITS + 1), None);
}

#[test]
fn wide_words_exceed_the_old_128_qubit_ceiling() {
    // The exact case that used to fail with "supports n_qubits ≤ 128".
    let (h, jumps) = model(130, 4);
    let spec = LindbladSpec::<W256>::new(130, &h, &jumps).unwrap();
    assert_eq!(spec.n_qubits(), 130);

    let (h, jumps) = model(512, 4);
    let spec = LindbladSpec::<W512>::new(512, &h, &jumps).unwrap();
    assert_eq!(spec.n_qubits(), 512);
}

#[test]
fn too_many_qubits_for_the_width_is_rejected() {
    let (h, jumps) = model(200, 4);
    let Err(err) = LindbladSpec::<W128>::new(200, &h, &jumps) else {
        panic!("a 200-qubit spec must not fit 128-qubit words");
    };
    assert_eq!(
        err.to_string(),
        "LindbladSpec supports n_qubits ≤ 128; got 200"
    );
    assert!(LindbladSpec::<W256>::new(200, &h, &jumps).is_ok());
    assert!(word_from_codes::<W128>(&[0u8; 129]).is_err());
}

#[test]
fn padding_into_a_wider_word_changes_nothing() {
    // Identical 6-qubit physics, embedded in 64-, 200- and 400-qubit
    // registers backed by 128-, 256- and 512-qubit words.
    let narrow = evolve::<W128>(64, 6, 6, false);
    let wide = evolve::<W256>(200, 6, 6, false);
    let widest = evolve::<W512>(400, 6, 6, false);
    assert!(narrow.abs() > 1e-6, "test observable is trivially zero");
    assert_eq!(narrow.to_bits(), wide.to_bits(), "{narrow} vs {wide}");
    assert_eq!(narrow.to_bits(), widest.to_bits(), "{narrow} vs {widest}");

    // The Kossakowski accumulation sums through a hash map, whose iteration
    // order follows the word's hash — and a wider word hashes differently.
    // Same physics, so agreement to rounding.
    let narrow = evolve::<W128>(64, 6, 6, true);
    for wide in [
        evolve::<W256>(200, 6, 6, true),
        evolve::<W512>(400, 6, 6, true),
    ] {
        assert!(
            (narrow - wide).abs() <= 1e-13 * narrow.abs(),
            "{narrow} vs {wide}"
        );
    }
}

#[test]
fn same_width_different_register_size_agrees() {
    // Within one width, the spectator qubits must not touch the answer.
    assert_eq!(
        evolve::<W256>(130, 6, 5, false).to_bits(),
        evolve::<W256>(256, 6, 5, false).to_bits()
    );
}

#[test]
fn orbit_rep_step_is_width_independent() {
    // The momentum-orbit path (canonicalization, masked-shift generators,
    // character table) on the same ring stored in 128- and 512-qubit words.
    // Untruncated, so both widths hold the same reps; basis *order* follows
    // hash-map iteration (width-dependent), so coefficients agree to rounding.
    for (n, k, steps) in [(12, 0, 3), (12, 1, 3), (100, 0, 2)] {
        let narrow = evolve_orbit::<W128>(n, k, steps);
        let wide = evolve_orbit::<W512>(n, k, steps);
        assert!(narrow.len() > 10, "basis did not grow");
        assert_eq!(narrow.len(), wide.len(), "n={n}, k={k}");
        for ((wa, ca), (wb, cb)) in narrow.iter().zip(&wide) {
            assert_eq!(wa, wb, "n={n}, k={k}: rep sets differ");
            assert!((ca - cb).norm() <= 1e-12, "n={n}, k={k}: {ca} vs {cb}");
        }
    }
}

#[test]
fn orbit_rep_step_runs_beyond_128_qubits() {
    let reps = evolve_orbit::<W256>(130, 0, 2);
    assert!(reps.len() > 10);
}
