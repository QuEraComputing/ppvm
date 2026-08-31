// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use fxhash::FxHashMap;
use num::Complex;

/// Test-only full-space complex predictor-corrector step, UNTRUNCATED:
/// adds every nonzero leakage string (two hops) and applies the exact
/// in-basis exponential to the complex coefficient vector. Reference
/// bridge between the real `pc_step` and the orbit-rep path.
fn pc_step_complex_full(
    spec: &LindbladSpec,
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    dt: f64,
) {
    let protected: Vec<Word> = Vec::new();
    let leak = spec.leakage_complex(basis, coeffs, &protected).unwrap();
    for (w, v) in leak {
        if v.norm() > 0.0 {
            basis.push(w);
            coeffs.push(Complex::new(0.0, 0.0));
        }
    }
    let predict = mf_expm::expm_apply_mf_cxvec(spec, basis, dt, coeffs, 0.0);
    let leak2 = spec.leakage_complex(basis, &predict, &protected).unwrap();
    drop(predict);
    for (w, v) in leak2 {
        if v.norm() > 0.0 {
            basis.push(w);
            coeffs.push(Complex::new(0.0, 0.0));
        }
    }
    *coeffs = mf_expm::expm_apply_mf_cxvec(spec, basis, dt, coeffs, 0.0);
}

fn jump_hpauli(s: &str, rate: f64) -> JumpInput {
    JumpInput {
        lincomb: vec![(s.to_string(), Complex::new(1.0, 0.0))],
        rate,
    }
}

#[test]
fn z_dephasing_action_on_x() {
    // L = Z on a single qubit; L*(X) = γ(ZXZ - X) = γ(-X - X) = -2γ X.
    let spec = LindbladSpec::new(
        1,
        &[("X".to_string(), 0.0)], // no Hamiltonian
        &[jump_hpauli("Z", 0.5)],
    )
    .unwrap();
    let (x, _) = parse_pauli_string("X", 1).unwrap();
    let terms = spec.action(&x);
    assert_eq!(terms.len(), 1);
    assert!((terms[0].1 - (-1.0)).abs() < 1e-12); // -2·0.5 = -1
}

#[test]
fn amplitude_damping_action_on_z() {
    // Single-qubit σ⁻ jump: L*(Z) = -γ(I + Z). With γ=1 we expect
    // I coefficient = -1, Z coefficient = -1.
    let sigma_minus = JumpInput {
        lincomb: vec![
            ("X".to_string(), Complex::new(0.5, 0.0)),
            ("Y".to_string(), Complex::new(0.0, -0.5)),
        ],
        rate: 1.0,
    };
    let spec = LindbladSpec::new(1, &[], &[sigma_minus]).unwrap();
    let (z, _) = parse_pauli_string("Z", 1).unwrap();
    let terms = spec.action(&z);
    let (i_word, _) = parse_pauli_string("I", 1).unwrap();
    let mut i_coeff = 0.0;
    let mut z_coeff = 0.0;
    for (w, c) in &terms {
        if w == &i_word {
            i_coeff = *c;
        } else if w == &z {
            z_coeff = *c;
        }
    }
    assert!((i_coeff - (-1.0)).abs() < 1e-10, "I coeff = {i_coeff}");
    assert!((z_coeff - (-1.0)).abs() < 1e-10, "Z coeff = {z_coeff}");
}

#[test]
fn word_codec_roundtrip() {
    let codes = [0u8, 1, 2, 3, 1, 0, 3, 2];
    let w = word_from_codes(&codes).unwrap();
    let mut out = vec![0u8; codes.len()];
    codes_from_word(&w, &mut out);
    assert_eq!(out.as_slice(), &codes);
}

/// The full-space complex step at momentum k=0 must reproduce the real
/// pc_step on the same trajectory exactly.
#[test]
fn complex_full_matches_real_at_kzero() {
    let n = 4usize;
    let dt = 0.01f64;
    let n_steps = 5usize;
    let mut h_terms: Vec<(String, f64)> = Vec::new();
    for j in 0..n {
        let nxt = (j + 1) % n;
        for op in ["X", "Y"] {
            let mut s = vec!['I'; n];
            s[j] = op.chars().next().unwrap();
            s[nxt] = op.chars().next().unwrap();
            h_terms.push((s.into_iter().collect(), 1.0));
        }
    }
    let spec = LindbladSpec::new(n, &h_terms, &[]).unwrap();

    let mut basis_r: Vec<Word> = (0..n)
        .map(|j| {
            let mut s = vec!['I'; n];
            s[j] = 'Z';
            let st: String = s.into_iter().collect();
            let (w, _) = parse_pauli_string(&st, n).unwrap();
            w
        })
        .collect();
    let mut coeffs_r: Vec<f64> = vec![1.0; n];

    let mut basis_c = basis_r.clone();
    let mut coeffs_c: Vec<Complex<f64>> = coeffs_r.iter().map(|&v| Complex::new(v, 0.0)).collect();

    let protected: Vec<Word> = Vec::new();
    for _ in 0..n_steps {
        // Large max_basis: rank cap never binds, so the real path
        // enriches fully (adds every leakage string). Match the
        // complex path by setting its tau_add=0.0 (also full
        // enrichment) so the two stay in lock-step at k=0.
        spec.pc_step(
            &mut basis_r,
            &mut coeffs_r,
            dt,
            &protected,
            &PcStepConfig {
                max_basis: 10_000_000,
                ..Default::default()
            },
        )
        .unwrap();
        pc_step_complex_full(&spec, &mut basis_c, &mut coeffs_c, dt);
    }
    // Match as (word → coeff) maps.
    let map_r: FxHashMap<Word, f64> = basis_r.into_iter().zip(coeffs_r).collect();
    let map_c: FxHashMap<Word, Complex<f64>> = basis_c.into_iter().zip(coeffs_c).collect();
    assert_eq!(
        map_r.len(),
        map_c.len(),
        "real and complex pc_step produced different basis sizes ({} vs {})",
        map_r.len(),
        map_c.len()
    );
    let mut max_diff = 0.0_f64;
    for (w, cr) in &map_r {
        let cc = map_c
            .get(w)
            .copied()
            .unwrap_or_else(|| panic!("word {:?} in real but not complex", w));
        assert!(cc.im.abs() < 1e-10, "expected zero imag at k=0, got {cc:?}");
        max_diff = max_diff.max((cr - cc.re).abs());
    }
    assert!(
        max_diff < 1e-10,
        "real vs complex pc_step diverged: max |Δc| = {max_diff:e}"
    );
}

/// Small-system end-to-end check that orbit-rep merging gives the
/// same physics as standard evolution, when no truncation is applied.
///
/// Setup: n=4 qubit chain, PBC, translation-invariant XY Hamiltonian
/// `H = Σ_j (X_j X_{j+1} + Y_j Y_{j+1})`, no dissipation. Initial
/// operator `O(0) = Σ_j Z_j` is translation-invariant (k=0 sector).
///
/// Run 10 pc_step iterations with `drop_tol = 0` (no truncation):
/// once without merging, once applying `canonicalize_pauli_sum`
/// after each step. Canonicalize the un-merged final state once at
/// the end. The two orbit-rep representations should be
/// bit-identical up to FP noise.
#[test]
fn pc_step_matches_symmetry_merged_on_small_chain() {
    use ppvm_pauli_sum::symmetry::{TranslationGroup, canonicalize_pauli_sum};

    let n = 4usize;
    let dt = 0.05f64;
    let n_steps = 10usize;

    // Build XY-chain Hamiltonian with PBC. 8 terms (4 bonds × {XX, YY}).
    let mut h_terms: Vec<(String, f64)> = Vec::new();
    for j in 0..n {
        let nxt = (j + 1) % n;
        for op in ["X", "Y"] {
            let mut s = vec!['I'; n];
            s[j] = op.chars().next().unwrap();
            s[nxt] = op.chars().next().unwrap();
            h_terms.push((s.into_iter().collect(), 1.0));
        }
    }
    // No dissipation.
    let spec = LindbladSpec::new(n, &h_terms, &[]).unwrap();
    let group = TranslationGroup::chain_1d(n);

    // Initial: O(0) = Σ_j Z_j (translation-invariant).
    let mut basis_u: Vec<Word> = (0..n)
        .map(|j| {
            let mut s = vec!['I'; n];
            s[j] = 'Z';
            let st: String = s.into_iter().collect();
            let (w, _) = parse_pauli_string(&st, n).unwrap();
            w
        })
        .collect();
    let mut coeffs_u: Vec<f64> = vec![1.0; n];

    // Mirror state for the "with merging" run.
    let mut basis_m = basis_u.clone();
    let mut coeffs_m = coeffs_u.clone();

    let protected: Vec<Word> = Vec::new();
    for _ in 0..n_steps {
        // max_basis == current basis size → room = 0: no leakage
        // enrichment, only the expm step (the regime where merging
        // commutes with evolution). drop_tol = 0 → no truncation.
        let cfg_u = PcStepConfig {
            max_basis: basis_u.len(),
            ..Default::default()
        };
        spec.pc_step(&mut basis_u, &mut coeffs_u, dt, &protected, &cfg_u)
            .unwrap();

        let cfg_m = PcStepConfig {
            max_basis: basis_m.len(),
            ..Default::default()
        };
        spec.pc_step(&mut basis_m, &mut coeffs_m, dt, &protected, &cfg_m)
            .unwrap();
        // Apply symmetry merging on the "with merging" run only.
        canonicalize_pauli_sum(&mut basis_m, &mut coeffs_m, &group);
    }

    // Canonicalize the un-merged final state once.
    canonicalize_pauli_sum(&mut basis_u, &mut coeffs_u, &group);

    // Both representations should now be in orbit-rep form; compare
    // as (word → coeff) maps with FP tolerance.
    let map_u: FxHashMap<Word, f64> = basis_u.into_iter().zip(coeffs_u).collect();
    let map_m: FxHashMap<Word, f64> = basis_m.into_iter().zip(coeffs_m).collect();
    assert_eq!(
        map_u.len(),
        map_m.len(),
        "merged basis size {} != post-merged-unmerged basis size {}",
        map_m.len(),
        map_u.len()
    );
    let mut max_diff = 0.0f64;
    for (w, c_u) in &map_u {
        let c_m = map_m.get(w).copied().unwrap_or_else(|| {
            panic!(
                "rep {:?} present in un-merged-then-canonicalized but not in merged",
                w
            );
        });
        max_diff = max_diff.max((c_u - c_m).abs());
    }
    assert!(
        max_diff < 1e-9,
        "with-merging vs without-merging diverged: max |Δc| = {max_diff:e}"
    );
}
