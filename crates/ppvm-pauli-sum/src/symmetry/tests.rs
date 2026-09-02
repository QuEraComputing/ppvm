// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use crate::sum::PauliSum;
use fxhash::FxHashMap;
use num::Complex;
use ppvm_pauli_word::word::PauliWord;
use std::f64::consts::PI;

type W = PauliWord<[u8; 1], fxhash::FxBuildHasher, true>;

fn word(s: &str) -> W {
    W::from(s)
}

#[test]
fn chain_1d_canonicalizes_via_cyclic_shift() {
    let g = TranslationGroup::chain_1d(4);
    // All cyclic shifts of "IIXY" should canonicalize to the same rep.
    let candidates = ["IIXY", "IXYI", "XYII", "YIIX"];
    let canon: Vec<W> = candidates
        .iter()
        .map(|s| g.canonicalize(&word(s)))
        .collect();
    for c in &canon[1..] {
        assert_eq!(
            *c, canon[0],
            "all cyclic shifts must canonicalize to same rep"
        );
    }
}

#[test]
fn chain_1d_canonicalize_is_lex_min() {
    let g = TranslationGroup::chain_1d(4);
    let canon = g.canonicalize(&word("YIIX"));
    let orbit: Vec<W> = g.orbit(&word("YIIX")).collect();
    let min = orbit.iter().min().unwrap();
    assert_eq!(canon, *min);
}

#[test]
fn orbit_has_correct_size_for_chain() {
    let g = TranslationGroup::chain_1d(4);
    // "XIII" has orbit of size 4 (full chain).
    let orbit: Vec<W> = g.orbit(&word("XIII")).collect();
    assert_eq!(orbit.len(), 4);
    // "XIXI" has orbit of size 2 (period-2 invariant); 4 elements
    // total in the orbit iterator, but only 2 unique.
    let orbit: Vec<W> = g.orbit(&word("XIXI")).collect();
    assert_eq!(orbit.len(), 4); // iterator yields |G|, including duplicates
    let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
    assert_eq!(unique.len(), 2);
}

#[test]
fn torus_2d_canonicalize() {
    // 3x2 torus, 6 qubits.
    let g = TranslationGroup::torus_2d(3, 2);
    assert_eq!(g.n_qubits(), 6);
    assert_eq!(g.order(), 6);
    // X at (0,0) — orbit is all 6 single-X positions.
    let w = word("XIIIII");
    let orbit: Vec<W> = g.orbit(&w).collect();
    let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
    assert_eq!(unique.len(), 6);
    // All canonicalize to the same rep.
    let canon = g.canonicalize(&w);
    for u in &unique {
        assert_eq!(g.canonicalize(u), canon);
    }
}

#[test]
fn ladder_canonicalize() {
    // 2-leg ladder, L=3 → 6 qubits, group order 3 (no swap of legs).
    let g = TranslationGroup::ladder(3, 2);
    assert_eq!(g.n_qubits(), 6);
    assert_eq!(g.order(), 3);
    // X on leg 0 site 0: orbit = {(0,0), (0,1), (0,2)}, NOT including leg 1 sites.
    let w = word("XIIIII"); // qubit 0 = X
    let orbit: Vec<W> = g.orbit(&w).collect();
    assert_eq!(orbit.len(), 3);
    let unique: std::collections::HashSet<W> = orbit.into_iter().collect();
    assert_eq!(unique.len(), 3);
    // The orbit should be {qubit 0=X, qubit 1=X, qubit 2=X} — all leg 0.
    let expected: std::collections::HashSet<W> = ["XIIIII", "IXIIII", "IIXIII"]
        .iter()
        .map(|s| word(s))
        .collect();
    assert_eq!(unique, expected);
}

#[test]
fn canonicalize_pauli_sum_merges_orbit_members() {
    let g = TranslationGroup::chain_1d(4);
    let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
    let mut coeffs: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
    canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
    // All four collapse to one rep with coeff 1+2+3+4 = 10.
    assert_eq!(basis.len(), 1);
    assert!((coeffs[0] - 10.0).abs() < 1e-12);
}

#[test]
fn canonicalize_pauli_sum_keeps_distinct_orbits() {
    let g = TranslationGroup::chain_1d(4);
    // Two distinct orbits: {XIII, ...} (size 4) and {ZIII, ...} (size 4).
    let mut basis: Vec<W> = vec![word("XIII"), word("IXII"), word("ZIII"), word("IZII")];
    let mut coeffs: Vec<f64> = vec![1.0, 1.0, 2.0, 2.0];
    canonicalize_pauli_sum(&mut basis, &mut coeffs, &g);
    assert_eq!(basis.len(), 2);
    // Coefficients should be {2.0, 4.0} in some order.
    let mut cs = coeffs.clone();
    cs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!((cs[0] - 2.0).abs() < 1e-12);
    assert!((cs[1] - 4.0).abs() < 1e-12);
}

#[test]
fn canonicalize_with_shift_round_trip() {
    // For each cyclic shift of "IIXY" by `a` positions, the shift
    // counter returned should reproduce the original word when
    // applied to the canonical rep.
    let g = TranslationGroup::chain_1d(4);
    for src in ["IIXY", "IXYI", "XYII", "YIIX"] {
        let w = word(src);
        let (rep, cnt) = g.canonicalize_with_shift(&w);
        // Apply gen 0 `cnt[0]` times to rep, should equal w.
        let mut cur = rep;
        for _ in 0..cnt[0] {
            cur = g.apply_generator(&cur, 0);
        }
        assert_eq!(cur, w, "shift {cnt:?} doesn't reproduce {src}");
    }
}

#[test]
fn canonicalize_in_sector_agrees_with_canonicalize_with_shift() {
    let g = TranslationGroup::chain_1d(4);
    for src in ["IIXY", "IXYI", "XYII", "YIIX", "XIXI", "IIII"] {
        let w = word(src);
        let (rep, shift, orbit_size) = g.canonicalize_in_sector(&w, &[0]).unwrap();
        let (ref_rep, ref_shift) = g.canonicalize_with_shift(&w);
        assert_eq!(rep, ref_rep, "{src}: rep");
        assert_eq!(shift, ref_shift, "{src}: shift");
        let distinct: std::collections::HashSet<W> = g.orbit(&w).collect();
        assert_eq!(orbit_size, distinct.len(), "{src}: orbit size");
    }
}

#[test]
fn canonicalize_in_sector_rejects_incompatible_stabilizer() {
    // "XIXI" has period 2 on a 4-site chain: 2 distinct orbit members,
    // stabilizer generated by T². χ_k(T²) = e^{iπk}, so the orbit
    // carries the k=0 and k=2 sectors but not k=1 or k=3.
    let g = TranslationGroup::chain_1d(4);
    let w = word("XIXI");
    for k in [0, 2] {
        let (_, _, orbit_size) = g
            .canonicalize_in_sector(&w, &[k])
            .unwrap_or_else(|| panic!("k={k} must be compatible with a period-2 orbit"));
        assert_eq!(orbit_size, 2, "k={k}");
    }
    for k in [1, 3] {
        assert!(
            g.canonicalize_in_sector(&w, &[k]).is_none(),
            "k={k} must be rejected on a period-2 orbit"
        );
    }
    // A free orbit carries every sector, with the full |G| members.
    for k in 0..4 {
        let (_, _, orbit_size) = g.canonicalize_in_sector(&word("XIII"), &[k]).unwrap();
        assert_eq!(orbit_size, 4, "k={k}");
    }
}

#[test]
fn character_trivial_sector_is_one() {
    let g = TranslationGroup::chain_1d(4);
    // k=0 mode → character is always 1.
    for cnt in [vec![0u32], vec![1u32], vec![2u32], vec![3u32]] {
        let chi = g.character(&[0], &cnt);
        assert!((chi - Complex::new(1.0, 0.0)).norm() < 1e-12);
    }
}

#[test]
fn character_obeys_unit_modulus() {
    let g = TranslationGroup::chain_1d(4);
    for k in 0..4 {
        for a in 0..4 {
            let chi = g.character(&[k], &[a as u32]);
            assert!(
                (chi.norm() - 1.0).abs() < 1e-12,
                "|χ_{k}(T^{a})| should be 1, got {}",
                chi.norm()
            );
        }
    }
}

#[test]
fn character_numerator_normalizes_negative_modes() {
    let group = TranslationGroup::chain_1d(4);
    assert_eq!(group.character_numerator(&[-1], &[1]), 3);
    assert_eq!(group.character_numerator(&[3], &[1]), 3);
    assert!((group.character(&[-1], &[1]) - Complex::new(0.0, -1.0)).norm() < 1e-12);
}

#[test]
fn exact_character_detects_cross_generator_kernel() {
    let swap = vec![1, 0];
    let group = TranslationGroup::from_generators(2, vec![swap.clone(), swap], vec![2, 2]);
    assert_ne!(group.character_numerator(&[1, 0], &[1, 1]), 0);
    assert_eq!(group.character_numerator(&[1, 1], &[1, 1]), 0);
}

#[test]
fn character_checks_slice_lengths_in_release_builds() {
    let group = TranslationGroup::chain_1d(4);
    assert!(std::panic::catch_unwind(|| group.character(&[], &[0])).is_err());
    assert!(std::panic::catch_unwind(|| group.character(&[0], &[])).is_err());
}

#[test]
fn period_two_k_zero_round_trip_preserves_rep_coefficient() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XIXI"), word("IXIX")];
    let mut coeffs = vec![Complex::new(1.0, 0.0); 2];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[0]);
    assert_eq!(basis.len(), 1);
    assert!((coeffs[0] - Complex::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn period_two_compatible_k_two_round_trip_preserves_rep_coefficient() {
    let group = TranslationGroup::chain_1d(4);
    let rep = group.canonicalize(&word("XIXI"));
    let mut members: FxHashMap<W, Complex<f64>> = FxHashMap::default();
    for (member, counter) in group.orbit_with_counters(&rep) {
        members
            .entry(member)
            .or_insert_with(|| group.character(&[2], &counter).conj());
    }
    let (mut basis, mut coeffs): (Vec<W>, Vec<Complex<f64>>) = members.into_iter().unzip();
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[2]);
    assert_eq!(basis, vec![rep]);
    assert!((coeffs[0] - Complex::new(1.0, 0.0)).norm() < 1e-12);
}

#[test]
fn incompatible_stabilizer_projects_orbit_to_zero() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XXXX")];
    let mut coeffs = vec![Complex::new(1.0, 0.0)];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[1]);
    assert!(basis.is_empty());
    assert!(coeffs.is_empty());
}

#[test]
fn partial_period_two_orbit_is_averaged_with_missing_member_zero() {
    let group = TranslationGroup::chain_1d(4);
    let mut basis = vec![word("XIXI")];
    let mut coeffs = vec![Complex::new(1.0, 0.0)];
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &group, &[0]);
    assert_eq!(basis.len(), 1);
    assert!((coeffs[0] - Complex::new(0.5, 0.0)).norm() < 1e-12);
}

#[test]
fn momentum_zero_complex_projection_is_orbit_average() {
    // k=0 sector: complex projection averages orbit members onto the rep;
    // plain canonicalize_pauli_sum sums all coefficients onto the rep.
    let g = TranslationGroup::chain_1d(4);
    let basis: Vec<W> = vec![word("XIII"), word("IXII"), word("IIXI"), word("IIIX")];
    let real_coeffs = vec![1.0, 2.0, 3.0, 4.0];

    let mut basis_real = basis.clone();
    let mut coeffs_real = real_coeffs.clone();
    canonicalize_pauli_sum(&mut basis_real, &mut coeffs_real, &g);

    let mut basis_c = basis.clone();
    let mut coeffs_c: Vec<Complex<f64>> =
        real_coeffs.iter().map(|&v| Complex::new(v, 0.0)).collect();
    canonicalize_pauli_sum_complex(&mut basis_c, &mut coeffs_c, &g, &[0]);

    // Plain merge sums all coefficients onto the single orbit-rep:
    // 1+2+3+4 = 10. Complex k=0 projection averages over the orbit
    // (size 4), so we expect 10/4 = 2.5 on the rep.
    assert_eq!(basis_real.len(), 1);
    assert_eq!(basis_c.len(), 1);
    assert!((coeffs_real[0] - 10.0).abs() < 1e-12);
    assert!((coeffs_c[0].re - 2.5).abs() < 1e-12);
    assert!(coeffs_c[0].im.abs() < 1e-12);
}

#[test]
fn sector_check_rejects_missing_orbit_members() {
    let group = TranslationGroup::chain_1d(4);
    let basis = vec![word("ZIII")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    assert!(matches!(
        check_momentum_sector(&basis, &coeffs, &group, &[0], 1e-12),
        Err(SectorCheckError::CoefficientMismatch { .. })
    ));
}

#[test]
fn sector_check_rejects_incompatible_stabilizer() {
    let group = TranslationGroup::chain_1d(4);
    let basis = vec![word("XXXX")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    assert!(matches!(
        check_momentum_sector(&basis, &coeffs, &group, &[1], 1e-12),
        Err(SectorCheckError::IncompatibleStabilizer { .. })
    ));
}

#[test]
fn sector_check_rejects_invalid_numeric_inputs() {
    let group = TranslationGroup::chain_1d(2);
    let basis = vec![word("ZI")];
    assert!(matches!(
        check_momentum_sector(&basis, &[Complex::new(1.0, 0.0)], &group, &[0], f64::NAN),
        Err(SectorCheckError::InvalidTolerance { .. })
    ));
    assert!(matches!(
        check_momentum_sector(&basis, &[Complex::new(f64::NAN, 0.0)], &group, &[0], 1e-12),
        Err(SectorCheckError::NonFiniteCoefficient { .. })
    ));
}

#[test]
fn sector_check_rejects_nonfinite_coalesced_coefficient() {
    let group = TranslationGroup::chain_1d(1);
    let basis = vec![word("X"), word("X")];
    let coeffs = vec![Complex::new(f64::MAX, 0.0); 2];
    assert!(matches!(
        check_momentum_sector(&basis, &coeffs, &group, &[0], 1e-12),
        Err(SectorCheckError::NonFiniteCoefficient { .. })
    ));
}

#[test]
fn sector_error_display_names_the_words() {
    let group = TranslationGroup::chain_1d(2);
    let basis = vec![word("ZI")];
    let coeffs = vec![Complex::new(1.0, 0.0)];
    let message = check_momentum_sector(&basis, &coeffs, &group, &[0], 1e-12)
        .unwrap_err()
        .to_string();
    assert!(message.contains("ZI") || message.contains("IZ"));
}

#[test]
fn momentum_eigenstate_check_passes() {
    // O = Σ_j e^{ikj} Z_j for k = 2π/4 (mode 1) is a momentum-k
    // eigenstate. check_momentum_sector should accept.
    let g = TranslationGroup::chain_1d(4);
    let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
    let k_mode: i32 = 1;
    // Sector condition: c_{T^a p} = e^{-2πi k a / N} c_p.
    // Picking c_{Z_0} = 1: c_{Z_a} = e^{-2πi · 1 · a / 4} = (-i)^a.
    let coeffs: Vec<Complex<f64>> = (0..4_i32)
        .map(|a| Complex::from_polar(1.0, -2.0 * PI * (k_mode as f64) * (a as f64) / 4.0))
        .collect();
    let res = check_momentum_sector(&basis, &coeffs, &g, &[k_mode], 1e-10);
    assert!(
        res.is_ok(),
        "valid k-eigenstate failed sector check: {res:?}"
    );
}

#[test]
fn momentum_eigenstate_check_fails_for_wrong_sector() {
    // Same eigenstate as above, but check against the wrong momentum.
    let g = TranslationGroup::chain_1d(4);
    let basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
    let coeffs: Vec<Complex<f64>> = (0..4_i32)
        .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
        .collect();
    // Check against k=0 (constant) — should fail.
    let res = check_momentum_sector(&basis, &coeffs, &g, &[0], 1e-10);
    assert!(res.is_err(), "k=1 eigenstate wrongly passed as k=0 sector");
}

#[test]
fn momentum_eigenstate_round_trip_merge_preserves_rep_coeff() {
    // Merge a k=1 eigenstate; the orbit-rep coefficient should be
    // unchanged (= 1.0 for our chosen normalization, picking
    // c_{Z_0} = 1).
    let g = TranslationGroup::chain_1d(4);
    let mut basis: Vec<W> = vec![word("ZIII"), word("IZII"), word("IIZI"), word("IIIZ")];
    let mut coeffs: Vec<Complex<f64>> = (0..4_i32)
        .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / 4.0))
        .collect();
    canonicalize_pauli_sum_complex(&mut basis, &mut coeffs, &g, &[1]);
    assert_eq!(basis.len(), 1);
    // The canonical rep of single-Z orbit is Z_0 (lex-min of
    // {ZIII, IZII, IIZI, IIIZ} is IIIZ since 'I' < 'Z' lex-wise on
    // the (xbits, zbits) tuple; let's just check we got a single
    // entry with norm 1.
    assert!(
        (coeffs[0].norm() - 1.0).abs() < 1e-10,
        "expected |c_rep|=1, got {}",
        coeffs[0].norm()
    );
}

/// Trotter-mode end-to-end check that `PauliSum::symmetry_merge`
/// matches plain Trotter evolution post-canonicalized.
///
/// Setup: n=4 qubit chain, PBC, XY rotations on each bond. Initial
/// operator `O(0) = Σ_j Z_j` is translation-invariant.
///
/// **dt must be tiny.** First-order Trotter on a chain with PBC is
/// only translation-equivariant up to `O(dt^2)` (gate-order
/// commutator errors are NOT themselves T-symmetric). The
/// "merge-after-each-step" trajectory and the "merge-at-end"
/// trajectory therefore diverge by an amount proportional to that
/// Trotter error. We test in the dt → 0 limit where the divergence
/// is below FP noise.
#[test]
fn pauli_sum_symmetry_merge_matches_plain_trotter() {
    use crate::config::indexmap::ByteFxHashF64;
    use crate::prelude::*;

    type Cfg = ByteFxHashF64<1>;

    let n: usize = 4;
    // Tiny dt — Trotter per-step error scales as dt^2 and shows up
    // as a translation-non-equivariant correction; we want it below
    // FP noise at the tolerance we assert below (1e-7).
    let dt = 1e-5_f64;
    let n_steps = 2usize;
    let group = TranslationGroup::chain_1d(n);

    // Total-Z initial: O(0) = Σ_j Z_j (translation-invariant).
    let mut o_u: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    let mut o_m: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    for j in 0..n {
        let mut s: Vec<char> = vec!['I'; n];
        s[j] = 'Z';
        let st: String = s.into_iter().collect();
        o_u += (st.as_str(), 1.0);
        o_m += (st.as_str(), 1.0);
    }
    assert_eq!(o_u.len(), n);
    assert_eq!(o_m.len(), n);

    // Apply XY Trotter steps to both copies. With merging, call
    // symmetry_merge_pauli_sum after each step.
    for _ in 0..n_steps {
        for j in 0..n {
            let nxt = (j + 1) % n;
            o_u.rxx(j, nxt, dt);
            o_u.ryy(j, nxt, dt);
            o_m.rxx(j, nxt, dt);
            o_m.ryy(j, nxt, dt);
        }
        symmetry_merge_pauli_sum(&mut o_m, &group);
    }

    // Canonicalize the un-merged result once at the end.
    symmetry_merge_pauli_sum(&mut o_u, &group);

    // Compare as (word → coeff) maps, FP tolerance.
    let u: FxHashMap<_, f64> = o_u.iter().map(|(w, c)| (*w, *c)).collect();
    let m: FxHashMap<_, f64> = o_m.iter().map(|(w, c)| (*w, *c)).collect();
    assert_eq!(
        u.len(),
        m.len(),
        "post-merge basis sizes differ: u={} vs m={}",
        u.len(),
        m.len()
    );
    let mut max_diff = 0.0_f64;
    for (w, &cu) in &u {
        let cm = *m.get(w).unwrap_or_else(|| {
            panic!("rep present in u but not in m: {:?}", w);
        });
        max_diff = max_diff.max((cu - cm).abs());
    }
    // At dt = 1e-5 over 2 steps, accumulated Trotter
    // commutator-induced T-eq error is ~2·dt^2·|H|^2 ≈ 1e-9; we
    // assert 1e-7 to leave safety margin.
    assert!(
        max_diff < 1e-7,
        "Trotter with-merging diverged from without-merging: max |Δc| = {max_diff:e}"
    );
}

/// Build the `(re, im)` real pair of the momentum-`k` eigenstate
/// `O_k = Σ_a e^{-2πi k a / n} Z_a` on an `n`-site chain.
fn seed_z_momentum_pair<Cfg>(n: usize, k: i32) -> (PauliSum<Cfg>, PauliSum<Cfg>)
where
    Cfg: ppvm_traits::Config<Coeff = f64>,
    PauliSum<Cfg>: for<'s> std::ops::AddAssign<(&'s str, f64)>,
{
    let mut re: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    let mut im: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    for a in 0..n {
        let mut s: Vec<char> = vec!['I'; n];
        s[a] = 'Z';
        let st: String = s.into_iter().collect();
        let phase = -2.0 * PI * (k as f64) * (a as f64) / (n as f64);
        re += (st.as_str(), phase.cos());
        im += (st.as_str(), phase.sin());
    }
    (re, im)
}

/// At `k = 0` the phase-aware pair merge must agree with the independent
/// real-coefficient `symmetry_merge_pauli_sum` code path — for *every*
/// orbit, free or stabilized. This is the
/// regression test for the summing-vs-averaging convention: a global
/// `|G|` rescale of the averaged projector agrees only on free orbits.
#[test]
fn momentum_merge_pair_matches_symmetry_merge_at_k0() {
    use crate::config::indexmap::ByteFxHashF64;
    use crate::prelude::*;

    type Cfg = ByteFxHashF64<1>;
    let n = 4usize;
    let group = TranslationGroup::chain_1d(n);

    let mut reference: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    let mut re: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    let mut im: PauliSum<Cfg> = PauliSum::builder().n_qubits(n).build();
    // Σ_j Z_j and Σ_j X_j X_{j+1}: free orbits (|orbit| = |G|).
    for j in 0..n {
        let mut z: Vec<char> = vec!['I'; n];
        z[j] = 'Z';
        let mut xx: Vec<char> = vec!['I'; n];
        xx[j] = 'X';
        xx[(j + 1) % n] = 'X';
        for (s, c) in [(z, 1.5), (xx, -0.25)] {
            let st: String = s.into_iter().collect();
            reference += (st.as_str(), c);
            re += (st.as_str(), c);
        }
    }
    // Orbits WITH a stabilizer, where a global |G| rescale would be wrong:
    // "ZZZZ" is translation-invariant (|orbit| = 1) and "ZIZI" has period 2.
    for (st, c) in [("ZZZZ", 0.75), ("ZIZI", 2.0), ("IZIZ", -0.5)] {
        reference += (st, c);
        re += (st, c);
    }
    // `im` needs an entry to exist; a zero coefficient must not survive.
    im += ("IIII", 0.0);

    symmetry_merge_pauli_sum(&mut reference, &group);
    momentum_merge_pauli_sum_pair(&mut re, &mut im, &group, &[0]);

    let expected: FxHashMap<_, f64> = reference.iter().map(|(w, c)| (*w, *c)).collect();
    let got: FxHashMap<_, f64> = re.iter().map(|(w, c)| (*w, *c)).collect();
    assert_eq!(expected.len(), got.len(), "basis sizes differ");
    for (w, &c) in &expected {
        let g = *got.get(w).unwrap_or_else(|| panic!("missing rep {w:?}"));
        assert!(
            (c - g).abs() < 1e-12,
            "rep {w:?}: symmetry_merge gave {c}, momentum_merge gave {g}"
        );
    }
    assert_eq!(im.len(), 0, "a purely real input must leave `im` empty");
}

/// A momentum-`k` eigenstate folds onto a single orbit rep whose
/// coefficient is the *summing* projector `Σ_{p ∈ orbit} χ_k(g_p) · c_p`,
/// computed here directly from the group API.
#[test]
fn momentum_merge_pair_matches_summing_projector() {
    use crate::config::indexmap::ByteFxHashF64;

    type Cfg = ByteFxHashF64<1>;
    let n = 4usize;
    let group = TranslationGroup::chain_1d(n);

    for k in 0..n as i32 {
        let (mut re, mut im) = seed_z_momentum_pair::<Cfg>(n, k);
        // Coefficient of each orbit member before merging.
        let before: FxHashMap<W, Complex<f64>> = (0..n)
            .map(|a| {
                let mut s: Vec<char> = vec!['I'; n];
                s[a] = 'Z';
                let phase = -2.0 * PI * (k as f64) * (a as f64) / (n as f64);
                (
                    word(&s.into_iter().collect::<String>()),
                    Complex::from_polar(1.0, phase),
                )
            })
            .collect();

        momentum_merge_pauli_sum_pair(&mut re, &mut im, &group, &[k]);

        let got_re: FxHashMap<_, f64> = re.iter().map(|(w, c)| (*w, *c)).collect();
        let got_im: FxHashMap<_, f64> = im.iter().map(|(w, c)| (*w, *c)).collect();
        assert!(
            !got_re.is_empty() || !got_im.is_empty(),
            "k={k}: merged away"
        );

        let rep = group.canonicalize(&word(&{
            let mut s: Vec<char> = vec!['I'; n];
            s[0] = 'Z';
            s.into_iter().collect::<String>()
        }));
        // Σ over the orbit of χ_k(g) · c_{g·rep}.
        let mut expected = Complex::new(0.0, 0.0);
        for (member, counter) in group.orbit_with_counters(&rep) {
            expected += group.character(&[k], &counter) * before[&member];
        }
        let got = Complex::new(
            got_re.get(&rep).copied().unwrap_or(0.0),
            got_im.get(&rep).copied().unwrap_or(0.0),
        );
        assert!(
            (got - expected).norm() < 1e-12,
            "k={k}: rep {rep:?} expected {expected:?}, got {got:?}"
        );
    }
}

#[test]
#[should_panic(expected = "k_modes length 2 != number of generators 1")]
fn momentum_merge_pair_rejects_wrong_momentum_length() {
    use crate::config::indexmap::ByteFxHashF64;

    type Cfg = ByteFxHashF64<1>;
    let group = TranslationGroup::chain_1d(4);
    let (mut re, mut im) = seed_z_momentum_pair::<Cfg>(4, 0);
    momentum_merge_pauli_sum_pair(&mut re, &mut im, &group, &[0, 0]);
}

#[test]
#[should_panic(expected = "PauliSum qubit count 4 != group qubit count 3")]
fn momentum_merge_pair_rejects_qubit_count_mismatch() {
    use crate::config::indexmap::ByteFxHashF64;

    type Cfg = ByteFxHashF64<1>;
    let group = TranslationGroup::chain_1d(3);
    let (mut re, mut im) = seed_z_momentum_pair::<Cfg>(4, 0);
    momentum_merge_pauli_sum_pair(&mut re, &mut im, &group, &[0]);
}

#[test]
#[should_panic(expected = "generator 0 order must be nonzero")]
fn rejects_zero_generator_order() {
    TranslationGroup::from_generators(2, vec![vec![1, 0]], vec![0]);
}

#[test]
#[should_panic(expected = "declared order 4 != exact permutation order 2")]
fn rejects_inflated_generator_order() {
    TranslationGroup::from_generators(2, vec![vec![1, 0]], vec![4]);
}

#[test]
#[should_panic(expected = "generators 0 and 1 do not commute")]
fn rejects_noncommuting_generators() {
    let swap_01 = vec![1, 0, 2];
    let swap_12 = vec![0, 2, 1];
    TranslationGroup::from_generators(3, vec![swap_01, swap_12], vec![2, 2]);
}

#[test]
fn try_from_generators_reports_every_precondition() {
    use super::GroupError;
    /// `(n_qubits, perms, orders, expected error)`
    type Case = (usize, Vec<Vec<u32>>, Vec<u32>, GroupError);
    let cases: Vec<Case> = vec![
        (
            2,
            vec![vec![1, 0]],
            vec![2, 2],
            GroupError::LengthMismatch {
                perms: 1,
                orders: 2,
            },
        ),
        (
            3,
            vec![vec![1, 0]],
            vec![2],
            GroupError::PermutationLength {
                generator: 0,
                len: 2,
                n_qubits: 3,
            },
        ),
        (
            2,
            vec![vec![1, 5]],
            vec![2],
            GroupError::TargetOutOfRange {
                generator: 0,
                target: 5,
                n_qubits: 2,
            },
        ),
        (
            2,
            vec![vec![1, 1]],
            vec![2],
            GroupError::DuplicateTarget {
                generator: 0,
                target: 1,
            },
        ),
        (
            2,
            vec![vec![1, 0]],
            vec![0],
            GroupError::ZeroOrder { generator: 0 },
        ),
        (
            2,
            vec![vec![1, 0]],
            vec![4],
            GroupError::OrderMismatch {
                generator: 0,
                declared: 4,
                exact: 2,
            },
        ),
        (
            3,
            vec![vec![1, 0, 2], vec![0, 2, 1]],
            vec![2, 2],
            GroupError::NonCommuting { left: 0, right: 1 },
        ),
    ];
    for (n_qubits, perms, orders, expected) in cases {
        let err = TranslationGroup::try_from_generators(n_qubits, perms, orders)
            .expect_err("must be rejected");
        assert_eq!(err, expected);
    }
    // Valid input still constructs, and matches the panicking constructor.
    let group = TranslationGroup::try_from_generators(4, vec![vec![1, 2, 3, 0]], vec![4]).unwrap();
    assert_eq!(group.order(), TranslationGroup::chain_1d(4).order());
}

#[test]
fn rejects_zero_lattice_dimensions() {
    assert!(std::panic::catch_unwind(|| TranslationGroup::chain_1d(0)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::torus_2d(0, 2)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::torus_3d(2, 0, 2)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::ladder(2, 0)).is_err());
}

#[test]
fn rejects_dimension_product_overflow_before_allocation() {
    assert!(std::panic::catch_unwind(|| TranslationGroup::torus_2d(usize::MAX, 2)).is_err());
    assert!(std::panic::catch_unwind(|| TranslationGroup::ladder(usize::MAX, 2)).is_err());
}

#[test]
#[cfg(target_pointer_width = "64")]
#[should_panic(expected = "site count")]
fn rejects_site_count_outside_u32_addressable_range() {
    super::group::validate_site_count(u32::MAX as usize + 2, "test");
}

#[test]
fn odometer_yields_expected_counter_order() {
    let group = TranslationGroup::torus_2d(2, 3);
    let counters: Vec<Vec<u32>> = group
        .orbit_with_counters(&word("XIIIII"))
        .map(|(_, counter)| counter)
        .collect();
    assert_eq!(
        counters,
        vec![
            vec![0, 0],
            vec![1, 0],
            vec![0, 1],
            vec![1, 1],
            vec![0, 2],
            vec![1, 2],
        ],
    );
}

#[test]
fn traversal_matches_brute_force_composition() {
    let group = TranslationGroup::torus_2d(2, 3);
    let source = word("XYZIII");
    for (candidate, counter) in group.orbit_with_counters(&source) {
        let mut brute = source;
        for (g, &count) in counter.iter().enumerate() {
            for _ in 0..count {
                brute = group.apply_generator(&brute, g);
            }
        }
        assert_eq!(candidate, brute);
    }
}

#[test]
fn public_word_width_checks_are_not_debug_only() {
    let group = TranslationGroup::chain_1d(4);
    let short = word("XI");
    assert!(std::panic::catch_unwind(|| group.canonicalize(&short)).is_err());
    assert!(std::panic::catch_unwind(|| group.canonicalize_with_shift(&short)).is_err());
    assert!(std::panic::catch_unwind(|| group.orbit(&short)).is_err());
}

#[test]
fn trivial_group_checks_word_width() {
    let group = TranslationGroup::from_generators(4, vec![], vec![]);
    let short = word("XI");
    assert!(std::panic::catch_unwind(|| group.canonicalize(&short)).is_err());
    assert!(std::panic::catch_unwind(|| group.canonicalize_with_shift(&short)).is_err());
}

#[test]
fn rejects_group_order_overflow() {
    let orders = if usize::BITS == 64 {
        vec![u32::MAX, u32::MAX, u32::MAX]
    } else {
        vec![u32::MAX, u32::MAX]
    };
    assert!(std::panic::catch_unwind(|| { super::group::checked_group_order(&orders) }).is_err());
}
