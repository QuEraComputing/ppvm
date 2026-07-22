// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
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
fn momentum_zero_complex_merge_matches_real_merge() {
    // k=0 sector: complex merge with all-real input should give
    // real-valued orbit-rep coefficients equal to the plain
    // canonicalize_pauli_sum result.
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
    // 1+2+3+4 = 10. Complex merge does the same with a 1/|G|
    // prefactor, so we expect 10/4 = 2.5 on the rep.
    assert_eq!(basis_real.len(), 1);
    assert_eq!(basis_c.len(), 1);
    assert!((coeffs_real[0] - 10.0).abs() < 1e-12);
    assert!((coeffs_c[0].re - 2.5).abs() < 1e-12);
    assert!(coeffs_c[0].im.abs() < 1e-12);
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
fn rejects_group_order_overflow() {
    let orders = if usize::BITS == 64 {
        vec![u32::MAX, u32::MAX, u32::MAX]
    } else {
        vec![u32::MAX, u32::MAX]
    };
    assert!(std::panic::catch_unwind(|| { super::group::checked_group_order(&orders) }).is_err());
}
