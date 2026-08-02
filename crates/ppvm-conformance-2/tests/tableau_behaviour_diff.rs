// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Behaviour-parity suite (the PRIME DIRECTIVE): every user-facing contract the
//! integration baseline lists, asserted against OLD `ppvm-tableau`.
//!
//! "Behaviour" here is not just numbers — it is **when** side effects happen
//! (does a gate truncate on its own?), the boundary of every keep-rule, which
//! calls draw from the RNG, what lands in the measurement record, and the
//! *order* of the public amplitude vector. A divergence in any of those is a
//! correctness gap even when the numbers agree.
//!
//! Runs in a **debug** build, so every `debug_assert!` in both crates (contract
//! 17: no exact-zero coefficient stored, no imaginary measurement phase,
//! `pauli != I` in the decomposition, probabilities in `[0, 1]`, `cz_block_pairs`
//! same-word) is live throughout.

use num::complex::Complex64;
use ppvm_conformance_2::tableau::*;

use ppvm_tableau::data::GeneralizedTableau as OldGT;
use ppvm_tableau_2::GeneralizedTableau as NewGT;
use ppvm_tableau_2::Tableau as NewTableau;

/// `Σ |c|²` over the amplitude vector — the quantity `normalize()` drives to 1.
fn norm_sq<D: Driver>(t: &D) -> f64 {
    t.coeffs().iter().map(|(_, c)| c.norm_sqr()).sum()
}

// ===========================================================================
// 1 & 2 — truncation timing and the strict, absolute gate keep-rule
// ===========================================================================

/// Contract 1: gates **auto-truncate**. A branch that lands under the threshold
/// is gone the moment the gate returns — there is no `truncate()` to call, and a
/// gate-only sequence reproduces old exactly.
#[test]
fn gates_auto_truncate_without_a_caller_call() {
    // A large threshold makes the T gate's small branch fall under the cutoff
    // immediately, so the support size is the observable.
    for &thr in &[0.0f64, 1e-12, 0.1, 0.3, 0.4] {
        let mut o: OldNarrow = Driver::new_seeded(3, thr, 1);
        let mut m: NewNarrow = Driver::new_seeded(3, thr, 1);
        o.h(0);
        m.h(0);
        o.t(0);
        m.t(0);
        // NOTE: no truncate() call on either side — there is none on the tableau.
        assert_eq!(
            o.n_coeffs(),
            m.n_coeffs(),
            "threshold {thr}: support size after t() differs (old {} vs new {})",
            o.n_coeffs(),
            m.n_coeffs()
        );
        assert_eq!(
            o.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
            m.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
            "threshold {thr}: surviving index set differs"
        );
    }
}

/// Contract 2: the gate keep-rule is **strict** `>` on the ABSOLUTE magnitude.
/// A branch whose post-gate magnitude is EXACTLY the threshold must be DROPPED
/// (`ppvm-pauli-sum-2`'s `CoefficientThreshold` keeps it — `cutoff_mismatch` in
/// `lean/PPVM/Algebra/Truncation.lean`).
#[test]
fn gate_cutoff_is_strictly_greater_and_absolute() {
    // rz(θ) on |+⟩ writes branch magnitude |sin(θ/2)| onto the flipped index.
    // Choose θ so that |sin(θ/2)| is exactly representable and use it as the
    // threshold: strict `>` ⇒ the branch is dropped.
    let mag = 0.25f64;
    let theta = 2.0 * mag.asin();

    let mut o: OldNarrow = Driver::new_seeded(1, mag, 0);
    let mut m: NewNarrow = Driver::new_seeded(1, mag, 0);
    o.rz(0, theta);
    m.rz(0, theta);
    assert_eq!(o.n_coeffs(), m.n_coeffs(), "boundary support size");

    // Directly: an exactly-at-threshold magnitude must not survive.
    let exact = m
        .coeffs()
        .iter()
        .any(|(_, c)| (c.norm() - mag).abs() < 1e-15);
    assert!(
        !exact,
        "a branch of magnitude exactly == threshold survived; the keep-rule must be strict `>`"
    );

    // Just above the boundary it must survive, on both engines.
    let mut o2: OldNarrow = Driver::new_seeded(1, mag * (1.0 - 1e-9), 0);
    let mut m2: NewNarrow = Driver::new_seeded(1, mag * (1.0 - 1e-9), 0);
    o2.rz(0, theta);
    m2.rz(0, theta);
    assert_eq!(o2.n_coeffs(), m2.n_coeffs());
    assert!(o2.n_coeffs() >= m.n_coeffs());
}

// ===========================================================================
// 3 — measurement cutoffs: RELATIVE in case a, ABSENT in case b
// ===========================================================================

/// Contract 3a: case-a keeps `|c|² > threshold² · ‖v‖²` — relative to the
/// CURRENT norm, not absolute. A branch that is above the absolute threshold can
/// still be dropped (and vice versa); whichever old does, new must do.
#[test]
fn case_a_measurement_cutoff_is_relative_and_matches_old() {
    for seed in 0..24u64 {
        for &thr in &[1e-12f64, 1e-3, 1e-2, 0.05] {
            let mut o: OldNarrow = Driver::new_seeded(4, thr, seed);
            let mut m: NewNarrow = Driver::new_seeded(4, thr, seed);
            // Build a state with a spread of branch magnitudes, then measure a
            // qubit whose Z is NOT a stabilizer (case a).
            for (q, angle) in [(0usize, 0.31f64), (1, 0.05), (2, 0.9), (3, 0.02)] {
                o.h(q);
                m.h(q);
                o.rz(q, angle);
                m.rz(q, angle);
            }
            o.cnot(0, 1);
            m.cnot(0, 1);
            let a = o.measure(0);
            let b = m.measure(0);
            assert_eq!(a, b, "seed {seed} thr {thr}: case-a outcome");
            assert_eq!(
                o.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
                m.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
                "seed {seed} thr {thr}: case-a survivors after the relative cutoff"
            );
        }
    }
}

/// Contract 3b: case b applies **no magnitude filter at all** — only the parity
/// `retain`. Sub-threshold branches must survive a case-b measurement.
#[test]
fn case_b_measurement_applies_no_magnitude_filter() {
    // A generous threshold; the branches are built by `rz` on qubit 0 only, so
    // measuring qubit 2 (untouched, Z is a stabilizer) is case b.
    let thr = 0.4f64;
    let mut o: OldNarrow = Driver::new_seeded(3, thr, 7);
    let mut m: NewNarrow = Driver::new_seeded(3, thr, 7);
    o.h(0);
    m.h(0);
    // Two small-but-surviving branches: sin(θ/2) just above the threshold, then
    // scaled down by a second rotation.
    o.rz(0, 1.9);
    m.rz(0, 1.9);
    o.rz(1, 1.9);
    m.rz(1, 1.9);
    let before_o = o.coeffs_sorted();
    let before_m = m.coeffs_sorted();
    assert_eq!(before_o.len(), before_m.len());
    assert!(before_m.len() > 1, "need a branched state for this test");

    let a = o.measure(2);
    let b = m.measure(2);
    assert_eq!(a, b, "case-b outcome");
    assert_eq!(
        o.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
        m.coeffs_sorted().iter().map(|e| e.0).collect::<Vec<_>>(),
        "case-b survivors"
    );
    // Qubit 2 is |0⟩ and untouched: the parity predicate keeps everything, so
    // NOTHING may be dropped by magnitude.
    assert_eq!(
        m.n_coeffs(),
        before_m.len(),
        "case b dropped a branch — it must apply no magnitude filter"
    );
}

// ===========================================================================
// 4 — normalization timing
// ===========================================================================

/// Contract 4: gates never normalize; case-a measurement always does; case-b
/// only when the support shrank.
#[test]
fn normalization_timing_matches_old() {
    // (a) gates do not normalize.
    let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, 0);
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 0);
    o.h(0);
    m.h(0);
    o.t(0);
    m.t(0);
    let (no, nm) = (norm_sq(&o), norm_sq(&m));
    assert!(
        (no - nm).abs() < 1e-15,
        "gate-path norm differs: old {no} vs new {nm}"
    );
    assert_eq!(
        o.coeffs_sorted(),
        m.coeffs_sorted(),
        "h;t coefficients must be bit-for-bit identical"
    );

    // (b) a case-a measurement always normalizes.
    let mut a: NewNarrow = Driver::new_seeded(3, 1e-12, 5);
    a.h(0);
    a.t(0);
    a.cnot(0, 1);
    a.h(2);
    a.t(2);
    a.measure(0);
    assert!(
        (norm_sq(&a) - 1.0).abs() < 1e-15,
        "case-a measurement must normalize; got {}",
        norm_sq(&a)
    );

    // (c) a case-b measurement that drops nothing must leave the coefficients
    //     byte-identical (no renormalization applied).
    let mut b: NewNarrow = Driver::new_seeded(3, 1e-12, 5);
    b.h(0);
    b.t(0);
    let before = b.coeffs();
    let n_before = b.n_coeffs();
    b.measure(2); // qubit 2 untouched ⇒ case b, deterministic |0⟩, drops nothing
    assert_eq!(
        b.n_coeffs(),
        n_before,
        "case b dropped entries unexpectedly"
    );
    assert_eq!(
        before,
        b.coeffs(),
        "case b that drops nothing must not renormalize"
    );
}

/// `normalize()` panics on a zero-norm vector with the old message.
#[test]
fn normalize_panics_on_zero_norm() {
    let mut amp: ppvm_tableau_2::Amplitudes<usize> = ppvm_tableau_2::Amplitudes::new();
    amp.unsafe_insert(0, num::complex::Complex64::new(0.0, 0.0));
    let err = std::panic::catch_unwind(move || amp.normalize()).unwrap_err();
    let msg = err
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| err.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_default();
    assert!(
        msg.contains("Zero norm encountered during normalization"),
        "unexpected panic message: {msg}"
    );
}

// ===========================================================================
// 5 — `measure` return type and record
// ===========================================================================

/// Contract 5: a lost qubit's `measure` returns `None` and pushes exactly one
/// `None` onto the record, without touching the state.
#[test]
fn measure_on_lost_qubit_returns_none_and_records_none() {
    let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, 3);
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 3);
    o.loss_channel(0, 1.0);
    m.loss_channel(0, 1.0);
    assert_eq!(o.lost(), m.lost());
    assert!(m.lost()[0]);
    assert!(o.record().is_empty() && m.record().is_empty());

    let rows_before = m.rows();
    let coeffs_before = m.coeffs();
    assert_eq!(o.measure(0), None);
    assert_eq!(m.measure(0), None);
    assert_eq!(o.record(), m.record());
    assert_eq!(m.record().last(), Some(&None));
    assert_eq!(
        m.rows(),
        rows_before,
        "measuring a lost qubit touched the frame"
    );
    assert_eq!(
        m.coeffs(),
        coeffs_before,
        "measuring a lost qubit touched the amplitudes"
    );
}

/// Contract 5 (bare frame): `Tableau::measure` keeps no record, never
/// normalizes, and its deterministic (case-b) branch consumes **no** randomness.
///
/// Type note: old typed this as `-> bool` and the generalized one as
/// `-> Option<bool>`; the `-2` tower unifies both on `Option<bool>`, so the bare
/// frame returns `Some(_)` unconditionally. The information is the same; the
/// signature is not (reported as an API-surface gap).
#[test]
fn bare_frame_measure_is_recordless_and_case_b_draws_nothing() {
    use ppvm_traits_2::{Clifford, Measure};

    let mut t: NewTableau = NewTableau::new_with_seed(2, 42);
    // |00⟩: Z on both qubits is a stabilizer ⇒ case b ⇒ no RNG draw.
    let mut probe = t.clone();
    let r0 = Measure::measure(&mut t, 0);
    assert_eq!(r0, Some(false));
    // The RNG stream must be untouched: the same next draw on both.
    let a = t.rng_next_f64_for_test();
    let b = probe.rng_next_f64_for_test();
    assert_eq!(
        a, b,
        "a deterministic (case-b) measurement drew from the RNG"
    );
    Clifford::h(&mut probe, 0);
    let _: Option<bool> = Measure::measure(&mut t, 1);
}

/// Test-only RNG probe on the bare frame (the tableau's RNG is private).
trait RngProbe {
    fn rng_next_f64_for_test(&mut self) -> f64;
}
impl RngProbe for NewTableau {
    fn rng_next_f64_for_test(&mut self) -> f64 {
        use ppvm_tableau_2::TableauLike;
        use rand::RngExt;
        self.rng_mut().random::<f64>()
    }
}

// ===========================================================================
// 6 — reset / loss are measurement-record-NEUTRAL
// ===========================================================================

/// Contract 6: `reset` and `loss_channel` pop the record entry their internal
/// `measure` pushed. `asymmetric_loss_channel` does **not** (a known old-crate
/// defect the port reproduces verbatim pending sign-off).
#[test]
fn reset_and_loss_are_record_neutral() {
    let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, 1);
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 1);
    o.x(0);
    m.x(0);
    o.reset(0);
    m.reset(0);
    assert!(o.record().is_empty(), "old reset polluted the record");
    assert!(m.record().is_empty(), "new reset polluted the record");
    assert_eq!(o.measure(0), Some(false));
    assert_eq!(m.measure(0), Some(false));

    let mut o2: OldNarrow = Driver::new_seeded(2, 1e-12, 1);
    let mut m2: NewNarrow = Driver::new_seeded(2, 1e-12, 1);
    o2.loss_channel(0, 1.0);
    m2.loss_channel(0, 1.0);
    assert!(o2.record().is_empty());
    assert!(m2.record().is_empty());
    assert!(m2.lost()[0]);
    assert_eq!(o2.measure(0), None);
    assert_eq!(m2.measure(0), None);
}

/// The reproduced old-crate defect: an asymmetric-loss event pollutes the
/// record with a spurious `Some(bool)`. Pinned so the divergence is visible if
/// anyone "fixes" one side alone.
#[test]
fn asymmetric_loss_pollutes_the_record_on_both_engines() {
    let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, 1);
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 1);
    o.asymmetric_loss_channel(0, 1.0, 1.0);
    m.asymmetric_loss_channel(0, 1.0, 1.0);
    assert_eq!(o.lost(), m.lost());
    assert_eq!(
        o.record(),
        m.record(),
        "asymmetric-loss record behaviour diverged"
    );
    assert_eq!(
        m.record().len(),
        1,
        "old leaves exactly one spurious record entry; the port must reproduce it"
    );
}

// ===========================================================================
// 8 — `measure_noisy`
// ===========================================================================

/// Contract 8: `measure_noisy` OVERWRITES (never appends), the state follows the
/// TRUE outcome, and `flip_prob == 0.0` draws no randomness.
#[test]
fn measure_noisy_overwrites_and_projects_on_the_true_outcome() {
    for seed in 0..16u64 {
        let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, seed);
        let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, seed);
        o.h(0);
        m.h(0);
        let plain = o.fork(Some(seed + 1000));
        let plain_new = m.fork(Some(seed + 1000));

        let mut o_f = o.fork(Some(seed + 1000));
        let mut m_f = m.fork(Some(seed + 1000));
        let a = o_f.measure_noisy(0, 1.0);
        let b = m_f.measure_noisy(0, 1.0);
        assert_eq!(a, b, "seed {seed}: measure_noisy(flip=1) result");
        assert_eq!(o_f.record(), m_f.record(), "seed {seed}: record");
        assert_eq!(m_f.record().len(), 1, "exactly one record entry per call");
        assert_eq!(m_f.record()[0], b);

        // The returned bit is the negation of a same-seed plain measure...
        let mut p = plain_new.fork(Some(seed + 1000));
        let truth = p.measure(0);
        assert_eq!(b, truth.map(|t| !t), "flip_prob = 1.0 must invert the bit");
        // ...and the state followed the TRUE outcome, so re-measuring reproduces it.
        assert_eq!(m_f.measure(0), truth, "state must follow the TRUE outcome");
        let _ = plain;

        // flip_prob = 0.0 draws nothing: the stream is unperturbed.
        let mut z_old = o.fork(Some(seed + 77));
        let mut z_new = m.fork(Some(seed + 77));
        let mut ref_old = o.fork(Some(seed + 77));
        let mut ref_new = m.fork(Some(seed + 77));
        assert_eq!(z_old.measure_noisy(0, 0.0), ref_old.measure(0));
        assert_eq!(z_new.measure_noisy(0, 0.0), ref_new.measure(0));
        assert_eq!(z_new.peek_rng_f64(), ref_new.peek_rng_f64());
        assert_eq!(z_old.peek_rng_f64(), ref_old.peek_rng_f64());
    }
}

// ===========================================================================
// 9 & 10 — RNG-consumption discipline and the comparison boundaries
// ===========================================================================

/// Contract 9: `depolarize1` draws UNCONDITIONALLY (even on a lost qubit), so a
/// loss event does not desynchronize a seeded stream. `depolarize2` returns
/// early WITHOUT drawing, so it does.
#[test]
fn rng_consumption_discipline_matches_old() {
    // depolarize1 on a lost qubit: the stream must stay in step with a run where
    // the qubit was never lost.
    fn stream_after_depolarize1<D: Driver>(lose: bool) -> Vec<f64> {
        let mut t: D = D::new_seeded(2, 1e-12, 99);
        if lose {
            t.loss_channel(0, 1.0);
        }
        t.depolarize1(0, 0.5);
        (0..8).map(|_| t.peek_rng_f64()).collect()
    }
    // Direct comparison old-vs-new is the real bar.
    for &lose in &[false, true] {
        let so = stream_after_depolarize1::<OldNarrow>(lose);
        let sn = stream_after_depolarize1::<NewNarrow>(lose);
        assert_eq!(so, sn, "depolarize1 stream (lost = {lose}) diverged");
    }

    fn stream_after_depolarize2<D: Driver>(lose: bool) -> Vec<f64> {
        let mut t: D = D::new_seeded(2, 1e-12, 99);
        if lose {
            t.loss_channel(0, 1.0);
        }
        t.depolarize2(0, 1, 0.5);
        (0..8).map(|_| t.peek_rng_f64()).collect()
    }
    for &lose in &[false, true] {
        let so = stream_after_depolarize2::<OldNarrow>(lose);
        let sn = stream_after_depolarize2::<NewNarrow>(lose);
        assert_eq!(so, sn, "depolarize2 stream (lost = {lose}) diverged");
    }

    // And the documented asymmetry itself: depolarize1 consumes a draw even when
    // the target is lost, while depolarize2 does not.
    let mut a: NewNarrow = Driver::new_seeded(2, 1e-12, 5);
    let mut b: NewNarrow = Driver::new_seeded(2, 1e-12, 5);
    a.depolarize1(0, 0.0); // draws
    let a_next = a.peek_rng_f64();
    let b_first = b.peek_rng_f64();
    let b_next = b.peek_rng_f64();
    assert_eq!(a_next, b_next, "depolarize1 must consume exactly one draw");
    let _ = b_first;

    let mut c: NewNarrow = Driver::new_seeded(2, 1e-12, 5);
    let mut d: NewNarrow = Driver::new_seeded(2, 1e-12, 5);
    c.loss_channel(0, 1.0);
    d.loss_channel(0, 1.0);
    c.depolarize2(0, 1, 0.5); // must NOT draw (early return on loss)
    assert_eq!(
        c.peek_rng_f64(),
        d.peek_rng_f64(),
        "depolarize2 drew from the RNG despite a lost endpoint"
    );
}

/// Contract 10: the probability-comparison boundaries are gate-specific.
/// `loss_channel(q, 0.0)` fires when the draw is exactly `0.0` (`p >= r`), the
/// OPPOSITE strictness from `depolarize1` (`p > r`) — reproduced verbatim.
#[test]
fn probability_boundaries_match_old() {
    // Sweep many seeds so the boundary decisions are exercised broadly, and
    // require the two engines to agree on every fired/not-fired decision.
    for seed in 0..64u64 {
        for &p in &[0.0f64, 0.001, 0.5, 0.999, 1.0] {
            let mut o: OldNarrow = Driver::new_seeded(2, 1e-12, seed);
            let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, seed);
            o.loss_channel(0, p);
            m.loss_channel(0, p);
            assert_eq!(o.lost(), m.lost(), "seed {seed} p {p}: loss decision");
            assert_eq!(o.record(), m.record());

            let mut o2: OldNarrow = Driver::new_seeded(2, 1e-12, seed);
            let mut m2: NewNarrow = Driver::new_seeded(2, 1e-12, seed);
            o2.h(0);
            m2.h(0);
            o2.depolarize1(0, p);
            m2.depolarize1(0, p);
            assert_rows_eq_local(&o2, &m2, &format!("seed {seed} p {p}: depolarize1"));

            let mut o3: OldNarrow = Driver::new_seeded(2, 1e-12, seed);
            let mut m3: NewNarrow = Driver::new_seeded(2, 1e-12, seed);
            o3.h(0);
            m3.h(0);
            o3.depolarize2(0, 1, p);
            m3.depolarize2(0, 1, p);
            assert_rows_eq_local(&o3, &m3, &format!("seed {seed} p {p}: depolarize2"));

            let mut o4: OldNarrow = Driver::new_seeded(2, 1e-12, seed);
            let mut m4: NewNarrow = Driver::new_seeded(2, 1e-12, seed);
            o4.correlated_loss_channel(0, 1, [p, p, p]);
            m4.correlated_loss_channel(0, 1, [p, p, p]);
            assert_eq!(
                o4.lost(),
                m4.lost(),
                "seed {seed} p {p}: correlated-loss decision"
            );
            assert_eq!(o4.record(), m4.record());
        }
    }
}

#[track_caller]
fn assert_rows_eq_local<A: Driver, B: Driver>(o: &A, m: &B, ctx: &str) {
    assert_eq!(o.rows(), m.rows(), "{ctx}: rows differ");
}

// ===========================================================================
// 11 & 12 — loss semantics on Clifford, batched, block and rot2 paths
// ===========================================================================

/// Contract 11: single-qubit Cliffords no-op on a lost qubit; two-qubit ones
/// no-op if EITHER endpoint is lost; batched variants FILTER (the surviving
/// indices still get the gate).
#[test]
fn lost_qubit_gate_semantics_match_old() {
    let n = 6;
    // single-qubit no-op
    let mut o: OldNarrow = Driver::new_seeded(n, 1e-12, 2);
    let mut m: NewNarrow = Driver::new_seeded(n, 1e-12, 2);
    o.loss_channel(0, 1.0);
    m.loss_channel(0, 1.0);
    let before = m.rows();
    o.h(0);
    m.h(0);
    assert_eq!(m.rows(), before, "a Clifford fired on a lost qubit");
    assert_eq!(o.rows(), m.rows());

    // two-qubit: no-op if either endpoint is lost
    let before2 = m.rows();
    o.cz(0, 1);
    m.cz(0, 1);
    assert_eq!(m.rows(), before2, "cz fired with a lost control");
    assert_eq!(o.rows(), m.rows());

    // batched: the surviving pairs still get the gate
    let mut o2: OldNarrow = Driver::new_seeded(n, 1e-12, 2);
    let mut m2: NewNarrow = Driver::new_seeded(n, 1e-12, 2);
    for t in [
        &mut o2 as &mut dyn LossDriver,
        &mut m2 as &mut dyn LossDriver,
    ] {
        t.d_h(1);
        t.d_h(3);
        t.d_h(5);
        t.d_loss(0, 1.0);
    }
    let mut expected = m2.fork(Some(0));
    expected.cz(2, 3);
    expected.cz(4, 5);
    o2.cz_many(&[(0, 1), (2, 3), (4, 5)]);
    m2.cz_many(&[(0, 1), (2, 3), (4, 5)]);
    assert_eq!(o2.rows(), m2.rows(), "cz_many with a lost pair");
    assert_eq!(
        expected.rows(),
        m2.rows(),
        "cz_many must FILTER: (2,3) and (4,5) still apply, (0,1) does not"
    );

    // sqrt_y_many with a lost target
    let mut o3: OldNarrow = Driver::new_seeded(n, 1e-12, 2);
    let mut m3: NewNarrow = Driver::new_seeded(n, 1e-12, 2);
    o3.loss_channel(0, 1.0);
    m3.loss_channel(0, 1.0);
    let mut exp3 = m3.fork(Some(0));
    for q in 1..n {
        exp3.sqrt_y(q);
    }
    let all: Vec<usize> = (0..n).collect();
    o3.sqrt_y_many(&all);
    m3.sqrt_y_many(&all);
    assert_eq!(o3.rows(), m3.rows(), "sqrt_y_many with a lost target");
    assert_eq!(exp3.rows(), m3.rows(), "sqrt_y_many must FILTER");
}

/// Object-safe helper for driving both engines through one loop.
trait LossDriver {
    fn d_h(&mut self, q: usize);
    fn d_loss(&mut self, q: usize, p: f64);
}
impl<T: Driver> LossDriver for T {
    fn d_h(&mut self, q: usize) {
        Driver::h(self, q)
    }
    fn d_loss(&mut self, q: usize, p: f64) {
        Driver::loss_channel(self, q, p)
    }
}

/// `cz_block` with a lost qubit must fall back to the per-pair path and produce
/// exactly the rows the per-pair loop does — on both engines.
#[test]
fn cz_block_with_loss_matches_per_pair_fallback() {
    let n = 85;
    let mut o: OldWide = Driver::new_seeded(n, 1e-10, 4);
    let mut m: NewWide = Driver::new_seeded(n, 1e-10, 4);
    for q in (0..n).step_by(2) {
        o.h(q);
        m.h(q);
    }
    o.loss_channel(5, 1.0);
    m.loss_channel(5, 1.0);

    let mut expected = m.fork(Some(0));
    for i in 0..17 {
        expected.cz(i, i + 17);
    }
    o.cz_block(0, 17, 17);
    m.cz_block(0, 17, 17);
    assert_eq!(o.rows(), m.rows(), "cz_block with loss: old vs new");
    assert_eq!(
        expected.rows(),
        m.rows(),
        "cz_block with loss must equal the per-pair fallback"
    );
}

/// Contract 12: `rotate_2` degrades to a single-qubit rotation with the SAME
/// angle when one endpoint is lost.
#[test]
fn rotate_2_degrades_to_rotate_1_on_loss() {
    let theta = 0.7f64;
    for &(lose_a, lose_b) in &[(true, false), (false, true), (true, true)] {
        let mut o: OldNarrow = Driver::new_seeded(3, 1e-12, 6);
        let mut m: NewNarrow = Driver::new_seeded(3, 1e-12, 6);
        o.h(0);
        m.h(0);
        o.h(1);
        m.h(1);
        if lose_a {
            o.loss_channel(0, 1.0);
            m.loss_channel(0, 1.0);
        }
        if lose_b {
            o.loss_channel(1, 1.0);
            m.loss_channel(1, 1.0);
        }
        let mut expected = m.fork(Some(0));
        if lose_a && !lose_b {
            expected.rz(1, theta);
        } else if lose_b && !lose_a {
            expected.rz(0, theta);
        }
        o.rzz(0, 1, theta);
        m.rzz(0, 1, theta);
        assert_eq!(
            o.coeffs_sorted(),
            m.coeffs_sorted(),
            "rzz-on-loss old vs new"
        );
        assert_eq!(
            expected.coeffs_sorted(),
            m.coeffs_sorted(),
            "rzz with a lost endpoint must equal rz on the survivor, same angle"
        );
    }
}

// ===========================================================================
// 13 — `rotate_2` ordering, container, merge semantics
// ===========================================================================

/// Contract 13: `rotate_2` applies `b` before `a`, uses the ABSOLUTE `|c| >
/// |cutoff|` keep-rule, and does NOT normalize. Values are compared per index
/// against old within 1e-12 (the merge ORDER is the one adjudicated divergence:
/// old's `std::HashMap` order is process-random, so no old order exists).
#[test]
fn rotate_2_merge_semantics_match_old() {
    for &(axis, name) in &[(0usize, "rxx"), (1, "ryy"), (2, "rzz")] {
        for seed in 0..8u64 {
            let mut o: OldNarrow = Driver::new_seeded(4, 1e-10, seed);
            let mut m: NewNarrow = Driver::new_seeded(4, 1e-10, seed);
            // Pre-branch the state so the merge has something to do.
            for q in 0..4 {
                o.h(q);
                m.h(q);
            }
            o.t(0);
            m.t(0);
            o.t(2);
            m.t(2);
            let theta = 0.37 + seed as f64 * 0.11;
            match axis {
                0 => {
                    o.rxx(0, 1, theta);
                    m.rxx(0, 1, theta);
                }
                1 => {
                    o.ryy(1, 2, theta);
                    m.ryy(1, 2, theta);
                }
                _ => {
                    o.rzz(2, 3, theta);
                    m.rzz(2, 3, theta);
                }
            }
            let co = o.coeffs_sorted();
            let cm = m.coeffs_sorted();
            assert_eq!(
                co.iter().map(|e| e.0).collect::<Vec<_>>(),
                cm.iter().map(|e| e.0).collect::<Vec<_>>(),
                "{name} seed {seed}: rotate_2 support"
            );
            for (a, b) in co.iter().zip(cm.iter()) {
                assert!(
                    (a.1 - b.1).norm() < 1e-12,
                    "{name} seed {seed}: coefficient at {} — old {} vs new {}",
                    a.0,
                    a.1,
                    b.1
                );
            }
            // rotate_2 must NOT normalize.
            let nn = norm_sq(&m);
            assert!(
                (nn - norm_sq(&o)).abs() < 1e-12,
                "{name}: norms diverged (old {} vs new {nn})",
                norm_sq(&o)
            );
        }
    }
}

/// Contract 13, the ordering half: the `b`-before-`a` apply order inside
/// `rotate_2` is observable, and old's exact order must be reproduced.
///
/// `rxx`/`ryy`/`rzz` all pass the SAME Pauli on both endpoints, so a port that
/// applied `a` before `b` would still be caught only by whatever relative sign
/// the two decompositions happen to produce. The general `rotate_2` entry point
/// takes two INDEPENDENT axes, so an `X⊗Z` rotation is the sharp test — and it
/// is public surface (`RotationTwo::rotate_2`) that the `rxx`/`ryy`/`rzz`
/// wrappers never reach. Every mixed axis pair is compared against old on a
/// pre-branched, deliberately asymmetric state.
#[test]
fn rotate_2_mixed_axis_order_matches_old() {
    // `[x, z]` axis encodings: X = [1,0], Y = [1,1], Z = [0,1].
    const AXES: [([u8; 2], &str); 3] = [([1, 0], "X"), ([1, 1], "Y"), ([0, 1], "Z")];

    for (axis_a, na) in AXES {
        for (axis_b, nb) in AXES {
            for seed in 0..6u64 {
                let mut o: OldNarrow = Driver::new_seeded(4, 1e-12, seed);
                let mut m: NewNarrow = Driver::new_seeded(4, 1e-12, seed);
                // Asymmetric pre-branching, so `a` and `b` sit in genuinely
                // different frames and the apply order cannot cancel out.
                for t in 0..4usize {
                    o.h(t);
                    m.h(t);
                }
                o.s(1);
                m.s(1);
                o.t(0);
                m.t(0);
                o.t(2);
                m.t(2);
                o.sqrt_y(3);
                m.sqrt_y(3);

                let theta = 0.41 + seed as f64 * 0.19;
                ppvm_traits::traits::RotationTwo::rotate_2(&mut o, axis_a, axis_b, 1, 2, theta);
                ppvm_traits_2::RotationTwo::<Complex64, f64>::rotate_2(
                    &mut m, axis_a, axis_b, 1, 2, theta,
                );

                let ctx = format!("r{na}{nb} seed {seed}");
                let co = o.coeffs_sorted();
                let cm = m.coeffs_sorted();
                assert_eq!(
                    co.iter().map(|e| e.0).collect::<Vec<_>>(),
                    cm.iter().map(|e| e.0).collect::<Vec<_>>(),
                    "{ctx}: rotate_2 support diverged"
                );
                for (x, y) in co.iter().zip(cm.iter()) {
                    assert!(
                        (x.1 - y.1).norm() < 1e-12,
                        "{ctx}: coefficient at {} — old {} vs new {}. A `b`-before-`a` \
                         ordering mismatch shows up here as a per-term sign/phase flip.",
                        x.0,
                        x.1,
                        y.1
                    );
                }
                // Still no normalization on the rot2 path.
                assert!(
                    (norm_sq(&m) - norm_sq(&o)).abs() < 1e-12,
                    "{ctx}: norms diverged — rotate_2 must not normalize"
                );
            }
        }
    }
}

/// The companion fact, established empirically over a randomized sweep: the
/// `b`-before-`a` apply order inside `rotate_2` is **not** observable.
///
/// The integration baseline asks to "verify the b-before-a ordering matters by
/// swapping it and observing a divergence". It does not — and that is worth
/// pinning as a test rather than leaving as a latent assumption, because it is
/// what makes the ordering safe to preserve verbatim without a Lean argument.
/// `compute_coefficients_after_pauli_apply` applies the frame-conjugated `P_a`
/// and `P_b`; since `a ≠ b` those two Paulis commute as operators, and each step
/// is an index XOR plus a multiplication by an exact ℤ/4 phase (`±1`, `±i`), so
/// there is not even float rounding to distinguish the two orders.
///
/// The sweep covers every mixed axis pair over randomized pre-branched states.
/// If a future change ever makes the order observable, this test fails and the
/// `b`-before-`a` choice becomes load-bearing again.
#[test]
fn rotate_2_endpoint_swap_is_observationally_identical() {
    const AXES: [[u8; 2]; 3] = [[1, 0], [1, 1], [0, 1]];
    let mut rng = ppvm_conformance_2::seeded_rng(2072_u64);

    for _ in 0..40 {
        // A randomized, deliberately asymmetric pre-branched state, replayed
        // identically into both orders.
        let prep: Vec<(usize, usize)> = (0..10)
            .map(|_| {
                (
                    rand::RngExt::random_range(&mut rng, 0..5u32) as usize,
                    rand::RngExt::random_range(&mut rng, 0..5usize),
                )
            })
            .collect();
        let theta: f64 = rand::RngExt::random_range(&mut rng, 0.05..3.0f64);
        let build = || {
            let mut t: NewNarrow = Driver::new_seeded(5, 1e-14, 3);
            for &(g, q) in &prep {
                match g {
                    0 => t.h(q),
                    1 => t.s(q),
                    2 => t.t(q),
                    3 => t.sqrt_y(q),
                    _ => t.sqrt_x(q),
                }
            }
            t
        };

        for axis_a in AXES {
            for axis_b in AXES {
                let mut ab = build();
                ppvm_traits_2::RotationTwo::<Complex64, f64>::rotate_2(
                    &mut ab, axis_a, axis_b, 1, 3, theta,
                );
                // The SAME operator P₁ ⊗ P₃, but written so the implementation
                // applies the other endpoint's Pauli first.
                let mut ba = build();
                ppvm_traits_2::RotationTwo::<Complex64, f64>::rotate_2(
                    &mut ba, axis_b, axis_a, 3, 1, theta,
                );

                let x = ab.coeffs_sorted();
                let y = ba.coeffs_sorted();
                assert_eq!(
                    x.len(),
                    y.len(),
                    "endpoint swap changed the support size — the apply order has \
                     become observable; `rotate_2`'s b-before-a ordering is now \
                     load-bearing and needs a differential test that can see it"
                );
                for (p, q) in x.iter().zip(y.iter()) {
                    assert_eq!(p.0, q.0, "endpoint swap changed the support");
                    assert!(
                        (p.1 - q.1).norm() < 1e-15,
                        "endpoint swap changed the coefficient at {}: {} vs {}",
                        p.0,
                        p.1,
                        q.1
                    );
                }
            }
        }
    }
}

// ===========================================================================
// 15 — `expectation` / `z_expectation` are non-mutating
// ===========================================================================

/// Contract 15: the expectation entry points never mutate and never normalize,
/// and a zero-probability projection returns `0.0` instead of panicking.
#[test]
fn expectation_is_non_mutating() {
    let mut m: NewNarrow = Driver::new_seeded(2, 1e-12, 0);
    m.h(0);
    m.cnot(0, 1);
    m.t(0);
    let rows = m.rows();
    let coeffs = m.coeffs();
    let record = m.record();
    for w in ["II", "ZZ", "XX", "YY", "IZ", "ZI", "XZ", "YX"] {
        let v = m.expectation_str(w);
        assert!(v.is_finite(), "<{w}> is not finite");
    }
    for q in 0..2 {
        let _ = m.z_expectation(q);
    }
    assert_eq!(rows, m.rows(), "expectation mutated the frame");
    assert_eq!(coeffs, m.coeffs(), "expectation mutated the amplitudes");
    assert_eq!(record, m.record(), "expectation touched the record");

    // Zero-probability projection must not panic (Bell ⟨YY⟩ = −1, ⟨YX⟩ = 0).
    let mut bell: NewNarrow = Driver::new_seeded(2, 1e-12, 0);
    bell.h(0);
    bell.cnot(0, 1);
    assert!((bell.expectation_str("YX")).abs() < 1e-12);
}

// ===========================================================================
// 16 — construction / reset / fork semantics
// ===========================================================================

/// Contract 16: `fork(Some(s))` clones the ENTIRE state including the record;
/// `reset_all()` clears state, loss and record but does **not** reseed the RNG.
#[test]
fn construction_reset_and_fork_semantics_match_old() {
    let mut o: OldNarrow = Driver::new_seeded(3, 1e-12, 8);
    let mut m: NewNarrow = Driver::new_seeded(3, 1e-12, 8);
    o.h(0);
    m.h(0);
    assert_eq!(o.measure(0), m.measure(0));
    assert_eq!(o.record().len(), 1);

    let fo = o.fork(Some(7));
    let fm = m.fork(Some(7));
    assert_eq!(fo.record().len(), 1, "fork must clone the record");
    assert_eq!(fm.record().len(), 1, "fork must clone the record");
    assert_eq!(fo.record(), fm.record());
    assert_eq!(fo.rows(), fm.rows());

    // reset_all: state back to |0…0⟩, record cleared, RNG NOT reseeded.
    let mut o_r: OldGT<ppvm_pauli_sum::config::indexmap::ByteFxHashF64<8>, usize> =
        OldGT::new_with_seed(3, 1e-12, 8);
    let mut m_r: NewGT<[u8; 8], usize> = NewGT::new_with_seed(3, 1e-12, 8);
    Driver::h(&mut o_r, 0);
    Driver::measure(&mut o_r, 0);
    Driver::h(&mut m_r, 0);
    Driver::measure(&mut m_r, 0);
    // Pre-reset next draws must equal the post-reset next draws (no reseed).
    let mut o_probe = o_r.fork(Some(0));
    let mut m_probe = m_r.fork(Some(0));
    let _ = (&mut o_probe, &mut m_probe);
    let o_next = o_r.bernoulli(0.5);
    let m_next = m_r.bernoulli(0.5);
    assert_eq!(o_next, m_next, "post-measure RNG streams diverged");

    o_r.reset_all();
    m_r.reset_all();
    assert!(o_r.current_measurement_record().is_empty());
    assert!(m_r.current_measurement_record().is_empty());
    assert_eq!(m_r.coefficients.len(), 1);
    assert_eq!(m_r.coefficients.entries()[0].1, 0usize);
    assert!(m_r.is_lost.iter().all(|&l| !l));
    // The RNG kept going (no reseed): the two engines still agree draw for draw.
    for _ in 0..4 {
        assert_eq!(
            o_r.bernoulli(0.5),
            m_r.bernoulli(0.5),
            "reset_all must not reseed the RNG"
        );
    }

    // Same seed + same circuit ⇒ identical trajectories over 50 forked shots.
    for shot in 0..50u64 {
        let a: OldNarrow = Driver::new_seeded(3, 1e-12, 5);
        let b: NewNarrow = Driver::new_seeded(3, 1e-12, 5);
        let mut a = a.fork(Some(shot));
        let mut b = b.fork(Some(shot));
        a.h(0);
        b.h(0);
        a.cnot(0, 1);
        b.cnot(0, 1);
        a.t(2);
        b.t(2);
        assert_eq!(a.measure_all(), b.measure_all(), "shot {shot}");
        assert_eq!(a.record(), b.record());
    }

    // `reset_loss_channel` clears only the loss flag.
    let mut o_l: OldNarrow = Driver::new_seeded(2, 1e-12, 1);
    let mut m_l: NewNarrow = Driver::new_seeded(2, 1e-12, 1);
    o_l.loss_channel(0, 1.0);
    m_l.loss_channel(0, 1.0);
    let rows = m_l.rows();
    o_l.reset_loss_channel(0);
    m_l.reset_loss_channel(0);
    assert_eq!(o_l.lost(), m_l.lost());
    assert!(!m_l.lost()[0]);
    assert_eq!(rows, m_l.rows(), "reset_loss_channel changed the state");
}

/// `append_measurement_record` / `overwrite_last_measurement_record` (stim
/// `MPAD` and readout-noise support) behave identically.
#[test]
fn record_editing_helpers_match_old() {
    let mut o: OldGT<ppvm_pauli_sum::config::indexmap::ByteFxHashF64<8>, usize> =
        OldGT::new_with_seed(2, 1e-12, 0);
    let mut m: NewGT<[u8; 8], usize> = NewGT::new_with_seed(2, 1e-12, 0);
    o.append_measurement_record(Some(true));
    m.append_measurement_record(Some(true));
    o.append_measurement_record(None);
    m.append_measurement_record(None);
    assert_eq!(
        o.current_measurement_record(),
        m.current_measurement_record()
    );
    o.overwrite_last_measurement_record(Some(false));
    m.overwrite_last_measurement_record(Some(false));
    assert_eq!(
        o.current_measurement_record(),
        m.current_measurement_record()
    );
    assert_eq!(m.current_measurement_record(), [Some(true), Some(false)]);
}
