// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The user-facing behavioural contracts of `ppvm-tableau-2`.
//!
//! Every test here pins a contract the port must reproduce from `ppvm-tableau`:
//! truncation timing and the keep-rule boundary, normalization timing, the
//! measurement-record discipline, the RNG-consumption discipline of the
//! channels, loss semantics, and the golden expectation values.

use num::complex::Complex64;
use ppvm_pauli_word_2::PauliWord;
use ppvm_tableau_2::prelude::*;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

type Tab = GeneralizedTableau<usize>;
type WideTab = GeneralizedTableau<u128>;

fn word(s: &str) -> PauliWord {
    s.into()
}

fn rng(seed: u64) -> SmallRng {
    SmallRng::seed_from_u64(seed)
}

fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() < tol,
        "expected {expected}, got {actual} (|Δ| = {})",
        (actual - expected).abs()
    );
}

fn norm_sq(tab: &Tab) -> f64 {
    tab.coefficients
        .iter()
        .map(|(c, _)| c.norm_sqr())
        .sum::<f64>()
}

fn bell() -> Tab {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    tab.h(0);
    tab.cnot(0, 1);
    tab
}

// ─── Golden expectation values (contract 15) ──────────────────────────────

#[test]
fn single_qubit_expectations() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-12);
    assert_close(tab.expectation(&word("Z")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("X")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("I")), 1.0, 1e-12);

    let mut plus: Tab = GeneralizedTableau::new(1, 1e-12);
    plus.h(0);
    assert_close(plus.expectation(&word("X")), 1.0, 1e-12);
    assert_close(plus.expectation(&word("Z")), 0.0, 1e-12);
}

#[test]
fn bell_state_pauli_expectations() {
    let mut tab = bell();
    assert_close(tab.expectation(&word("II")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("ZZ")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("XX")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("YY")), -1.0, 1e-12);
    assert_close(tab.expectation(&word("IZ")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("ZI")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("XZ")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("YX")), 0.0, 1e-12);
}

#[test]
fn ghz_state_expectations() {
    let mut tab: Tab = GeneralizedTableau::new(3, 1e-12);
    tab.h(0);
    tab.cnot(0, 1);
    tab.cnot(1, 2);
    assert_close(tab.expectation(&word("III")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("ZZZ")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("ZIZ")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("ZZI")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("IZI")), 0.0, 1e-12);
    assert_close(tab.expectation(&word("XXX")), 1.0, 1e-12);
    assert_close(tab.expectation(&word("YYY")), 0.0, 1e-12);
}

#[test]
fn ry_rotation_expectations_track_the_angle() {
    for theta in [0.0, 0.3, 1.0, std::f64::consts::FRAC_PI_2] {
        let mut tab: Tab = GeneralizedTableau::new(1, 1e-12);
        tab.ry(0, theta);
        assert_close(tab.expectation(&word("Z")), theta.cos(), 1e-12);
        assert_close(tab.expectation(&word("X")), theta.sin(), 1e-12);
    }
}

/// After `H; T` the state is `(|0⟩ + e^{iπ/4}|1⟩)/√2`. `⟨Y⟩` in particular
/// drives `phase_decomp` odd, exercising the `phase == 1 | 3` arms of the
/// case-a overlap — a sign bug there would flip `⟨Y⟩`.
#[test]
fn t_gate_superposition_expectations() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-12);
    tab.h(0);
    tab.t(0);
    assert_close(
        tab.expectation(&word("X")),
        std::f64::consts::FRAC_1_SQRT_2,
        1e-12,
    );
    assert_close(
        tab.expectation(&word("Y")),
        std::f64::consts::FRAC_1_SQRT_2,
        1e-12,
    );
}

/// `expectation` never mutates: the amplitude vector and the frame are
/// byte-identical afterwards.
#[test]
fn expectation_does_not_mutate() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    tab.h(0);
    tab.t(0);
    tab.cnot(0, 1);
    let before_coeffs: Vec<_> = tab.coefficients.iter().copied().collect();
    let before_rows: Vec<_> = tab.tableau.rows().collect();

    for w in ["II", "ZZ", "XY", "YZ"] {
        let _ = tab.expectation(&word(w));
        let _ = tab.z_expectation(0);
    }

    assert_eq!(
        before_coeffs,
        tab.coefficients.iter().copied().collect::<Vec<_>>()
    );
    assert_eq!(before_rows, tab.tableau.rows().collect::<Vec<_>>());
}

// ─── Truncation timing and boundary (contracts 1, 2, 4) ───────────────────

/// Gates auto-truncate: no `truncate()` call is needed, and none exists.
/// The keep-rule is **strict** `>` on the absolute magnitude, so a branch whose
/// magnitude is *exactly* the threshold is dropped.
#[test]
fn gate_cutoff_is_strict_and_fires_inline() {
    // With threshold 0 nothing is dropped; read the small branch's magnitude.
    let mut open: Tab = GeneralizedTableau::new(1, 0.0);
    open.rx(0, 0.05);
    assert_eq!(open.coefficients.len(), 2, "rx must branch");
    let small = open
        .coefficients
        .iter()
        .map(|(c, _)| c.norm())
        .fold(f64::INFINITY, f64::min);
    assert!(small > 0.0);

    // Threshold set exactly at that magnitude: strict `>` drops it.
    let mut at: Tab = GeneralizedTableau::new(1, small);
    at.rx(0, 0.05);
    assert_eq!(
        at.coefficients.len(),
        1,
        "a branch exactly at the threshold must be DROPPED (strict >)"
    );

    // Just below: kept.
    let mut below: Tab = GeneralizedTableau::new(1, small * (1.0 - 1e-9));
    below.rx(0, 0.05);
    assert_eq!(below.coefficients.len(), 2);
}

/// `rotate_2` uses a *third* rule: absolute on the magnitude itself
/// (`|c| > |threshold|`), still strict.
#[test]
fn rotate_2_cutoff_is_strict_absolute_magnitude() {
    let mut open: Tab = GeneralizedTableau::new(2, 0.0);
    open.rxx(0, 1, 0.05);
    assert_eq!(open.coefficients.len(), 2);
    let small = open
        .coefficients
        .iter()
        .map(|(c, _)| c.norm())
        .fold(f64::INFINITY, f64::min);

    let mut at: Tab = GeneralizedTableau::new(2, small);
    at.rxx(0, 1, 0.05);
    assert_eq!(at.coefficients.len(), 1);
}

/// Gates never normalize — not even the one that just dropped a branch. A
/// case-a measurement always does.
#[test]
fn normalization_timing() {
    // The threshold drops the small branch, so the norm leaves 1 and stays there.
    let mut tab: Tab = GeneralizedTableau::new(1, 0.1);
    tab.rx(0, 0.05);
    assert_eq!(
        tab.coefficients.len(),
        1,
        "the small branch must be dropped"
    );
    let after_gate = norm_sq(&tab);
    assert!(
        (after_gate - 1.0).abs() > 1e-6,
        "a truncating gate must NOT renormalize (norm² = {after_gate})"
    );
    tab.rz(0, 0.3);
    assert_close(norm_sq(&tab), after_gate, 1e-15);

    let mut m: Tab = GeneralizedTableau::new(1, 1e-12);
    let mut m_rng = rng(3);
    m.h(0);
    m.t(0);
    let _ = m.measure(0, &mut m_rng); // case a: Z is not a stabilizer of H|0⟩
    assert_close(norm_sq(&m), 1.0, 1e-15);
}

/// A case-b measurement that drops nothing must leave the amplitudes
/// byte-identical — no renormalization, and the pre-existing order preserved.
#[test]
fn case_b_measurement_that_drops_nothing_is_byte_identical() {
    let mut tab: Tab = GeneralizedTableau::new(2, 0.1);
    let mut rng = rng(1);
    tab.h(0);
    tab.cnot(0, 1);
    tab.rx(0, 0.05); // truncates → norm != 1, and Z_0 stays non-deterministic
    let before: Vec<_> = tab.coefficients.iter().copied().collect();

    // Qubit 1 is perfectly correlated with 0; measuring 0 first makes 1
    // deterministic (case b) with nothing to drop.
    let r0 = tab.measure(0, &mut rng).unwrap();
    let mid: Vec<_> = tab.coefficients.iter().copied().collect();
    let r1 = tab.measure(1, &mut rng).unwrap();
    assert_eq!(r0, r1, "Bell outcomes must agree");
    assert_eq!(
        mid,
        tab.coefficients.iter().copied().collect::<Vec<_>>(),
        "a case-b measurement that drops nothing must not touch the amplitudes"
    );
    assert!(!before.is_empty());
}

/// The amplitude vector is in ascending index order after a branching gate —
/// the sort-merge output. The order is public behaviour.
#[test]
fn branching_gate_leaves_ascending_order() {
    let mut tab: Tab = GeneralizedTableau::new(4, 1e-12);
    for q in 0..4 {
        tab.h(q);
        tab.t(q);
    }
    let idx: Vec<usize> = tab.coefficients.iter().map(|(_, i)| *i).collect();
    assert!(idx.len() > 8, "the T layer must branch");
    assert!(
        idx.windows(2).all(|w| w[0] < w[1]),
        "amplitudes must be ascending by index, got {idx:?}"
    );
}

// ─── Measurement record discipline (contracts 5, 6, 7, 8) ─────────────────

#[test]
fn measure_returns_none_and_records_none_for_a_lost_qubit() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut rng = rng(1);
    tab.loss_channel(0, 1.0, &mut rng);
    assert!(tab.is_lost[0]);
    assert!(
        tab.current_measurement_record().is_empty(),
        "a loss event is measurement-record-neutral"
    );

    assert_eq!(tab.measure(0, &mut rng), None);
    assert_eq!(tab.current_measurement_record(), &[None]);
}

#[test]
fn reset_is_measurement_record_neutral() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-12);
    let mut rng = rng(0);
    tab.x(0);
    tab.reset(0, &mut rng);
    assert!(tab.current_measurement_record().is_empty());
    assert_eq!(tab.measure(0, &mut rng), Some(false));
}

/// `measure_all`, `measure_many(all)` and a per-qubit `measure` loop must be
/// observationally identical — outcomes *and* record *and* RNG-draw order.
#[test]
fn batched_measurement_matches_a_per_qubit_loop() {
    let mut base: Tab = GeneralizedTableau::new(6, 1e-10);
    let mut setup_rng = rng(17);
    for q in 0..6 {
        base.h(q);
    }
    base.cnot(0, 1);
    base.t(2);
    base.cz(3, 4);
    base.t(4);
    base.loss_channel(5, 1.0, &mut setup_rng); // one lost qubit, so `None` entries appear

    let mut a = base.fork();
    let mut b = base.fork();
    let mut c = base.fork();
    let mut ar = rng(7);
    let mut br = rng(7);
    let mut cr = rng(7);

    let ra = a.measure_all(&mut ar);
    let rb = b.measure_many(&(0..6).collect::<Vec<_>>(), &mut br);
    let rc: Vec<Option<bool>> = (0..6).map(|q| c.measure(q, &mut cr)).collect();

    assert_eq!(ra, rb);
    assert_eq!(ra, rc);
    assert!(ra.contains(&None), "the lost qubit must report None");
    assert_eq!(
        a.current_measurement_record(),
        b.current_measurement_record()
    );
    assert_eq!(
        a.current_measurement_record(),
        c.current_measurement_record()
    );
}

/// A caller-owned scratch threaded across shots must not change any outcome.
#[test]
fn scratch_reuse_across_shots_is_observationally_neutral() {
    let mut base: Tab = GeneralizedTableau::new(5, 1e-10);
    for q in 0..5 {
        base.h(q);
    }
    base.t(0);
    base.t(3);

    let mut scratch: MeasureScratch<usize> = MeasureScratch::new();
    let mut with: Vec<Vec<Option<bool>>> = Vec::new();
    let mut without: Vec<Vec<Option<bool>>> = Vec::new();
    for seed in 0..8 {
        let mut with_tab = base.fork();
        let mut without_tab = base.fork();
        let mut with_rng = rng(seed);
        let mut without_rng = rng(seed);
        with.push(with_tab.measure_all_with_scratch(&mut scratch, &mut with_rng));
        without.push(without_tab.measure_all(&mut without_rng));
    }
    assert_eq!(with, without);
}

#[test]
fn measure_noisy_overwrites_the_record_and_keeps_the_true_projection() {
    let mut plain: Tab = GeneralizedTableau::new(1, 1e-12);
    let mut plain_rng = rng(11);
    plain.h(0);
    let truth = plain.measure(0, &mut plain_rng).unwrap();

    let mut noisy: Tab = GeneralizedTableau::new(1, 1e-12);
    let mut noisy_rng = rng(11);
    noisy.h(0);
    let reported = noisy.measure_noisy(0, 1.0, &mut noisy_rng).unwrap();

    assert_eq!(reported, !truth, "flip_prob = 1 must invert the report");
    assert_eq!(
        noisy.current_measurement_record(),
        &[Some(reported)],
        "exactly one record entry, holding the reported bit"
    );
    // The state followed the TRUE outcome, so re-measuring reproduces it.
    assert_eq!(noisy.measure(0, &mut noisy_rng), Some(truth));
}

/// `flip_with_prob(bit, 0.0)` returns `bit` and draws nothing.
#[test]
fn zero_flip_probability_does_not_perturb_the_stream() {
    let mut a = rng(5);
    let mut b = rng(5);
    assert!(!Tab::flip_with_prob(false, 0.0, &mut a));
    let sa: Vec<bool> = (0..16).map(|_| a.random::<f64>() < 0.5).collect();
    let sb: Vec<bool> = (0..16).map(|_| b.random::<f64>() < 0.5).collect();
    assert_eq!(sa, sb);
}

#[test]
fn deterministic_measurement_and_flip_do_not_perturb_the_stream() {
    let mut measured: Tab = GeneralizedTableau::new(1, 1e-12);
    let mut actual = rng(6);
    let mut reference = rng(6);

    assert_eq!(measured.measure(0, &mut actual), Some(false));
    assert!(Tab::flip_with_prob(false, 1.0, &mut actual));
    assert_eq!(
        actual.random::<u64>(),
        reference.random::<u64>(),
        "deterministic operations consumed an RNG draw"
    );
}

// ─── RNG-consumption discipline (contracts 9, 10) ─────────────────────────

/// `depolarize1` draws nothing on a lost qubit.
#[test]
fn single_qubit_depolarize_skips_the_draw_on_loss() {
    let mut lost: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut lost_rng = rng(21);
    lost.is_lost[0] = true;
    lost.depolarize1(0, 0.5, &mut lost_rng);

    let mut present: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut present_rng = rng(21);
    present.depolarize1(0, 0.5, &mut present_rng);

    let sa: Vec<bool> = (0..24).map(|_| lost_rng.random::<f64>() < 0.5).collect();
    let sb: Vec<bool> = (0..24).map(|_| present_rng.random::<f64>() < 0.5).collect();
    assert_ne!(sa, sb, "the lost-target channel must consume no draw");
}

/// The two-qubit channels return early **without** drawing, so a loss event
/// *does* shift the stream. Reproduced from the old crate verbatim (see the
/// crate's `# Deferrals`).
#[test]
fn two_qubit_depolarize_skips_the_draw_on_loss() {
    let mut lost: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut lost_rng = rng(21);
    lost.is_lost[0] = true;
    lost.depolarize2(0, 1, 0.5, &mut lost_rng);

    let mut present: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut present_rng = rng(21);
    present.depolarize2(0, 1, 0.5, &mut present_rng);

    let sa: Vec<bool> = (0..24).map(|_| lost_rng.random::<f64>() < 0.5).collect();
    let sb: Vec<bool> = (0..24).map(|_| present_rng.random::<f64>() < 0.5).collect();
    assert_ne!(sa, sb);
}

/// `depolarize1(q, 0.0)` is a state and stream no-op.
#[test]
fn zero_probability_channels_are_state_and_stream_noops() {
    let mut a: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut a_rng = rng(2);
    let before: Vec<_> = a.tableau.rows().collect();
    a.depolarize1(0, 0.0, &mut a_rng);
    a.pauli_error(0, [0.0; 3], &mut a_rng);
    a.depolarize2(0, 1, 0.0, &mut a_rng);
    a.two_qubit_pauli_error(0, 1, [0.0; 15], &mut a_rng);
    a.loss_channel(0, 0.0, &mut a_rng);
    a.correlated_loss_channel(0, 1, [0.0; 3], &mut a_rng);
    assert_eq!(before, a.tableau.rows().collect::<Vec<_>>());
    assert!(!a.is_lost.iter().any(|lost| *lost));

    let mut b_rng = rng(2);
    let sa: Vec<bool> = (0..8).map(|_| a_rng.random::<f64>() < 0.5).collect();
    let sb: Vec<bool> = (0..8).map(|_| b_rng.random::<f64>() < 0.5).collect();
    assert_eq!(sa, sb, "p = 0 must consume no draw");
}

// ─── Loss semantics (contracts 11, 12) ────────────────────────────────────

#[test]
fn clifford_gates_no_op_on_lost_qubits() {
    let mut tab: Tab = GeneralizedTableau::new(4, 1e-12);
    tab.is_lost[1] = true;
    let before: Vec<_> = tab.tableau.rows().collect();

    tab.h(1);
    tab.s(1);
    tab.cnot(1, 2); // control lost
    tab.cz(2, 1); // target lost
    tab.cy(1, 3);
    assert_eq!(before, tab.tableau.rows().collect::<Vec<_>>());

    // A pair with neither endpoint lost still applies.
    tab.cnot(2, 3);
    assert_ne!(before, tab.tableau.rows().collect::<Vec<_>>());
}

/// A batched gate **filters** rather than skipping wholesale: the surviving
/// indices/pairs still get the gate.
#[test]
fn batched_gates_filter_lost_qubits() {
    let mut batched: Tab = GeneralizedTableau::new(6, 1e-12);
    batched.is_lost[0] = true;
    let mut manual = batched.clone();

    batched.cz_many(&[(0, 1), (2, 3), (4, 5)]);
    manual.cz(2, 3);
    manual.cz(4, 5);
    assert_eq!(
        batched.tableau.rows().collect::<Vec<_>>(),
        manual.tableau.rows().collect::<Vec<_>>()
    );

    let mut b2: Tab = GeneralizedTableau::new(6, 1e-12);
    b2.is_lost[3] = true;
    let mut m2 = b2.clone();
    b2.sqrt_y_many(&[1, 3, 5]);
    m2.sqrt_y(1);
    m2.sqrt_y(5);
    assert_eq!(
        b2.tableau.rows().collect::<Vec<_>>(),
        m2.tableau.rows().collect::<Vec<_>>()
    );
}

/// `rotate_2` with one lost endpoint degrades to the single-qubit rotation on
/// the survivor, with the **same** angle.
#[test]
fn rotate_2_degrades_to_rotate_1_on_loss() {
    let theta = 0.37;

    let mut two: Tab = GeneralizedTableau::new(2, 1e-12);
    two.h(1);
    two.is_lost[0] = true;
    two.rzz(0, 1, theta);

    let mut one: Tab = GeneralizedTableau::new(2, 1e-12);
    one.h(1);
    one.is_lost[0] = true;
    one.rz(1, theta);

    let a: Vec<_> = two.coefficients.iter().copied().collect();
    let b: Vec<_> = one.coefficients.iter().copied().collect();
    assert_eq!(a.len(), b.len());
    for ((ca, ia), (cb, ib)) in a.iter().zip(b.iter()) {
        assert_eq!(ia, ib);
        assert!((ca - cb).norm() < 1e-12);
    }
}

#[test]
fn reset_loss_channel_only_clears_the_flag() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut rng = rng(0);
    tab.loss_channel(0, 1.0, &mut rng);
    assert!(tab.is_lost[0]);
    let before: Vec<_> = tab.tableau.rows().collect();

    tab.reset_loss_channel(0);
    assert!(!tab.is_lost[0]);
    assert_eq!(before, tab.tableau.rows().collect::<Vec<_>>());
}

// ─── Batched Clifford equivalence ─────────────────────────────────────────

/// Every fused batch sweep must agree with the per-qubit loop it replaces.
#[test]
fn fused_batches_match_the_single_qubit_loops() {
    let idx = [0usize, 3, 7, 12, 40, 63, 70, 84];
    let pairs = [(0usize, 17usize), (5, 40), (30, 70), (63, 64)];

    macro_rules! check_single {
        ($many:ident, $one:ident) => {{
            let mut batched: Tab = GeneralizedTableau::new(64, 1e-12);
            let mut manual: Tab = GeneralizedTableau::new(64, 1e-12);
            for q in 0..64 {
                batched.h(q);
                manual.h(q);
                batched.t(q % 7);
                manual.t(q % 7);
            }
            let small: Vec<usize> = idx.iter().copied().filter(|&i| i < 64).collect();
            batched.$many(&small);
            for &q in &small {
                manual.$one(q);
            }
            assert_eq!(
                batched.tableau.rows().collect::<Vec<_>>(),
                manual.tableau.rows().collect::<Vec<_>>(),
                concat!(stringify!($many), " must match the loop")
            );
        }};
    }

    check_single!(x_many, x);
    check_single!(y_many, y);
    check_single!(z_many, z);
    check_single!(h_many, h);
    check_single!(s_many, s);
    check_single!(s_dag_many, s_dag);
    check_single!(sqrt_x_many, sqrt_x);
    check_single!(sqrt_x_dag_many, sqrt_x_dag);
    check_single!(sqrt_y_many, sqrt_y);
    check_single!(sqrt_y_dag_many, sqrt_y_dag);

    macro_rules! check_pair {
        ($many:ident, $one:ident) => {{
            let mut batched: WideTab = GeneralizedTableau::new(85, 1e-12);
            let mut manual: WideTab = GeneralizedTableau::new(85, 1e-12);
            for q in 0..85 {
                batched.h(q);
                manual.h(q);
            }
            batched.$many(&pairs);
            for &(c, t) in &pairs {
                manual.$one(c, t);
            }
            assert_eq!(
                batched.tableau.rows().collect::<Vec<_>>(),
                manual.tableau.rows().collect::<Vec<_>>(),
                concat!(stringify!($many), " must match the loop")
            );
        }};
    }

    check_pair!(cnot_many, cnot);
    check_pair!(cz_many, cz);
    check_pair!(cy_many, cy);
}

// ─── Frame-only tableau ───────────────────────────────────────────────────

/// A deterministic (case-b) frame measurement consumes no randomness: a
/// following random (case-a) measurement must give the same bit as if the
/// case-b measurement had not happened at all.
#[test]
fn frame_case_b_measurement_leaves_the_rng_untouched() {
    for seed in 0..32 {
        let mut a: Tableau = Tableau::new(2);
        let mut ar = rng(seed);
        assert_eq!(a.measure(0, &mut ar), Some(false)); // |0⟩: Z is a stabilizer → case b
        a.h(1);

        let mut b: Tableau = Tableau::new(2);
        let mut br = rng(seed);
        b.h(1);

        assert_eq!(
            a.measure(1, &mut ar),
            b.measure(1, &mut br),
            "case b must not draw (seed {seed})"
        );
    }
}

/// The frame's `reset` is measure-then-`X`; it leaves `|0⟩`.
#[test]
fn frame_reset_returns_to_zero() {
    let mut tab: Tableau = Tableau::new(2);
    let mut rng = rng(0);
    tab.h(0);
    tab.reset(0, &mut rng);
    assert_eq!(tab.measure(0, &mut rng), Some(false));
}

// ─── Construction, fork, determinism ──────────────────────────────────────

#[test]
fn fork_clones_the_record_and_logical_state() {
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut rng = rng(3);
    tab.h(0);
    let _ = tab.measure(0, &mut rng);
    assert_eq!(tab.current_measurement_record().len(), 1);

    let forked = tab.fork();
    assert_eq!(forked.current_measurement_record().len(), 1);
    assert_eq!(
        forked.tableau.rows().collect::<Vec<_>>(),
        tab.tableau.rows().collect::<Vec<_>>()
    );
}

/// Same seed + same circuit ⇒ identical trajectories over many forked shots.
#[test]
fn seeded_runs_are_reproducible() {
    let build = || {
        let mut t: Tab = GeneralizedTableau::new(4, 1e-10);
        t.h(0);
        t.cnot(0, 1);
        t.t(2);
        t.h(2);
        t.rxx(2, 3, 0.4);
        t
    };
    let a = build();
    let b = build();
    for shot in 0..50 {
        let mut ar = rng(shot);
        let mut br = rng(shot);
        assert_eq!(a.fork().measure_all(&mut ar), b.fork().measure_all(&mut br));
    }
}

/// A fresh state is `|0…0⟩` with a single unit amplitude and no loss.
#[test]
fn fresh_state_is_the_computational_zero() {
    let tab: Tab = GeneralizedTableau::new(3, 1e-12);
    assert_eq!(tab.coefficients.len(), 1);
    assert_eq!(tab.coefficients.get(&0), Complex64::new(1.0, 0.0));
    assert!(tab.is_lost.iter().all(|&l| !l));
    assert!(tab.current_measurement_record().is_empty());
}

// ─── Physics smoke tests (ported from the old crate) ──────────────────────

#[test]
fn rotations_on_the_computational_basis() {
    let mut rng = rng(0);
    let mut rx: Tab = GeneralizedTableau::new(1, 1e-12);
    rx.rx(0, std::f64::consts::PI);
    assert_eq!(rx.coefficients.len(), 1, "rx(π) must not branch");
    assert_eq!(rx.measure(0, &mut rng), Some(true));

    let mut rz: Tab = GeneralizedTableau::new(1, 1e-12);
    rz.rz(0, 0.123);
    assert_eq!(rz.coefficients.len(), 1);
    assert_eq!(rz.measure(0, &mut rng), Some(false));

    let mut ry: Tab = GeneralizedTableau::new(1, 1e-10);
    ry.ry(0, 2.0 * std::f64::consts::PI);
    assert_eq!(ry.coefficients.len(), 1);
    assert_eq!(ry.measure(0, &mut rng), Some(false));
}

#[test]
fn two_qubit_rotations_correlate_the_outcomes() {
    // rxx(π/2)|00⟩ = (|00⟩ − i|11⟩)/√2: the outcomes must agree.
    let mut tab: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut tab_rng = rng(0);
    tab.rxx(0, 1, std::f64::consts::FRAC_PI_2);
    assert_eq!(tab.coefficients.len(), 2);
    assert_eq!(tab.measure(0, &mut tab_rng), tab.measure(1, &mut tab_rng));

    // rxx(π/2)|10⟩ = (|10⟩ − i|01⟩)/√2: the outcomes must differ.
    let mut anti: Tab = GeneralizedTableau::new(2, 1e-12);
    let mut anti_rng = rng(0);
    anti.x(0);
    anti.rxx(0, 1, std::f64::consts::FRAC_PI_2);
    assert_eq!(anti.coefficients.len(), 2);
    assert_ne!(
        anti.measure(0, &mut anti_rng),
        anti.measure(1, &mut anti_rng)
    );

    // rzz never branches on a computational basis state.
    let mut rzz: Tab = GeneralizedTableau::new(2, 1e-12);
    rzz.rzz(0, 1, std::f64::consts::FRAC_PI_2);
    assert_eq!(rzz.coefficients.len(), 1);
}

#[test]
fn u3_matches_its_rz_ry_rz_decomposition() {
    let (theta, phi, lambda) = (0.34, 0.21, 0.46);
    let mut u3: Tab = GeneralizedTableau::new(1, 1e-12);
    u3.u3(0, theta, phi, lambda);

    let mut manual: Tab = GeneralizedTableau::new(1, 1e-12);
    manual.rz(0, lambda);
    manual.ry(0, theta);
    manual.rz(0, phi);

    for seed in 0..64 {
        let mut ur = rng(seed);
        let mut mr = rng(seed);
        assert_eq!(
            u3.fork().measure(0, &mut ur),
            manual.fork().measure(0, &mut mr)
        );
    }
}

#[test]
fn r_matches_rz_rx_rz() {
    let (axis_angle, theta) = (0.21, 0.34);
    let mut r: Tab = GeneralizedTableau::new(1, 1e-12);
    r.r(0, axis_angle, theta);

    let mut manual: Tab = GeneralizedTableau::new(1, 1e-12);
    manual.rz(0, -axis_angle);
    manual.rx(0, theta);
    manual.rz(0, axis_angle);

    for seed in 0..64 {
        let mut rr = rng(seed);
        let mut mr = rng(seed);
        assert_eq!(
            r.fork().measure(0, &mut rr),
            manual.fork().measure(0, &mut mr)
        );
    }
}

#[test]
fn t_followed_by_t_dag_is_the_identity() {
    let mut tab: Tab = GeneralizedTableau::new(1, 1e-12);
    tab.h(0);
    let before: Vec<_> = tab.coefficients.iter().copied().collect();
    tab.t(0);
    tab.t_dag(0);
    let after: Vec<_> = tab.coefficients.iter().copied().collect();
    assert_eq!(before.len(), after.len());
    for ((ca, ia), (cb, ib)) in before.iter().zip(after.iter()) {
        assert_eq!(ia, ib);
        assert!((ca - cb).norm() < 1e-12, "T·T† must be the identity");
    }
}

// ─── Wide (multi-word) smoke test ─────────────────────────────────────────

/// An MSD-shaped 85-qubit circuit: fused batches, a cross-word `cz_block`, a
/// branching T layer, and a full `measure_all` sweep. Pins the cross-word bit
/// addressing and the wide (`u128`) bitstring path.
#[test]
fn wide_msd_shaped_circuit_runs_and_agrees_with_the_naive_form() {
    fn build(fused: bool) -> WideTab {
        let mut tab: WideTab = GeneralizedTableau::new(85, 1e-10);
        let block0: Vec<usize> = (0..17).collect();
        let block1: Vec<usize> = (17..34).collect();
        if fused {
            tab.sqrt_y_many(&block0);
            tab.sqrt_x_many(&block1);
            tab.cz_block(0, 17, 17);
            tab.cz_block(34, 51, 17);
        } else {
            for &q in &block0 {
                tab.sqrt_y(q);
            }
            for &q in &block1 {
                tab.sqrt_x(q);
            }
            for i in 0..17 {
                tab.cz(i, 17 + i);
            }
            for i in 0..17 {
                tab.cz(34 + i, 51 + i);
            }
        }
        for q in 0..8 {
            tab.h(q);
            tab.t(q);
        }
        tab
    }

    let fused = build(true);
    let naive = build(false);
    assert_eq!(
        fused.tableau.rows().collect::<Vec<_>>(),
        naive.tableau.rows().collect::<Vec<_>>(),
        "the fused batch surface must leave an identical frame"
    );
    assert_eq!(fused.coefficients.len(), naive.coefficients.len());
    assert!(fused.coefficients.len() > 1, "the T layer must branch");

    for shot in 0..8 {
        let mut fused_rng = rng(shot);
        let mut naive_rng = rng(shot);
        assert_eq!(
            fused.fork().measure_all(&mut fused_rng),
            naive.fork().measure_all(&mut naive_rng)
        );
    }
}
