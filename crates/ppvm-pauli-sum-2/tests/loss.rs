// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The neutral-atom loss workload, ported term for term from
//! `ppvm-pauli-sum/tests/loss.rs`.
//!
//! Old compares with `PartialEq`, i.e. **exact map equality including
//! zero-coefficient entries**, which makes `test_reset_channel` the sharpest test
//! of the no-implicit-reduce contract: `reset_loss_channel` on an already-lost
//! site scales the coefficient to exactly `0.0` and the term must stay.

use ppvm_pauli_sum_2::{
    LossyPauliSum, LossyPauliWord, MaxLossWeight, NoPolicy, PauliPattern, SiteSet, Sum,
};
use ppvm_traits_2::{
    Clifford, CorrelatedLossChannel, LossChannel, ResetLossChannel, RotationOne, Trace,
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// The `-2` sum backends are density-matrix-like: the noise and loss channels
/// are analytic coefficient scalings and never draw. Call sites thread this
/// fixed-seed RNG only to satisfy the injected-RNG trait surface.
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(0)
}

type L1 = LossyPauliSum<f64, NoPolicy>;

fn lw(s: &str) -> LossyPauliWord {
    LossyPauliWord::from(s)
}

fn sum(n: usize, terms: &[(&str, f64)]) -> L1 {
    Sum::from_terms(n, terms.iter().map(|(k, v)| (lw(k), *v)))
}

fn zero_state() -> PauliPattern {
    PauliPattern::repeated(SiteSet::Z.union(SiteSet::I))
}

// --- reset_loss_channel (old `test_reset_channel`) --------------------------

#[test]
fn reset_channel_maps_i_and_z_to_a_lost_branch() {
    // X / Y: no-op.
    for p in ["X", "Y"] {
        let mut state = sum(1, &[(p, 1.0)]);
        let before = state.clone();
        state.reset_loss_channel(0);
        assert_eq!(state, before, "reset_loss_channel must not touch {p}");
    }

    // I / Z: the term stays and an `L` branch is added at the same coefficient.
    for p in ["I", "Z"] {
        let mut state = sum(1, &[(p, 1.0)]);
        let mut want = state.clone();
        state.reset_loss_channel(0);
        want += (lw("L"), 1.0);
        assert_eq!(state, want, "reset_loss_channel on {p}");
    }
}

/// The no-implicit-reduce anchor: an already-lost term is scaled to **exactly
/// zero and kept**, so the result equals `state.clone() *= 0.0` under exact-map
/// equality (old's `state2 *= 0.0`).
#[test]
fn reset_channel_zeroes_a_lost_term_but_keeps_it() {
    let mut state = sum(1, &[("L", 1.0)]);
    let mut want = state.clone();
    state.reset_loss_channel(0);
    want *= 0.0;
    assert_eq!(state, want);
    assert_eq!(state.len(), 1, "the zeroed term must stay in the support");
    assert!(state.contains(&lw("L"), &0.0));
}

/// The `L` branch **accumulates** onto a colliding target rather than replacing
/// it (old's `test_reset_loss_channel_accumulates_duplicate_target_*`).
#[test]
fn reset_channel_accumulates_duplicate_targets() {
    let mut state = sum(1, &[("I", 2.0), ("Z", 3.0)]);
    state.reset_loss_channel(0);
    assert!(state.contains(&lw("I"), &2.0));
    assert!(state.contains(&lw("Z"), &3.0));
    assert!(state.contains(&lw("L"), &5.0), "2.0 + 3.0, not 3.0");
}

// --- loss_channel (old `test_loss_channel`) ---------------------------------

#[test]
fn loss_channel_scales_present_sites_and_branches_lost_ones() {
    for p in ["X", "Y", "I", "Z"] {
        let mut state = sum(1, &[(p, 1.0)]);
        let mut want = state.clone();
        state.loss_channel(0, 0.2, &mut rng());
        want *= 0.8;
        assert_eq!(state, want, "loss_channel on {p}");
    }

    // A lost site branches to `I` at `p` and is itself left **unscaled**.
    let mut state = sum(1, &[("L", 1.0)]);
    let mut want = state.clone();
    state.loss_channel(0, 0.2, &mut rng());
    want += (lw("I"), 0.2);
    assert_eq!(state, want);
}

// --- end-to-end: the loss-interleaved circuits ------------------------------

/// Old's `test_single_qubit_loss`: reset → identity Cliffords → loss → `X`, then
/// contract against the zero-state pattern. The `L` term is **not** counted by
/// the pattern, so the scalar is `−0.9 + 0.1`.
#[test]
fn single_qubit_loss_overlap() {
    let mut state = sum(1, &[("Z", 1.0)]);
    state.reset_loss_channel(0);

    let intermediate = state.clone();
    state.x(0);
    state.x(0);
    assert_eq!(state, intermediate, "x∘x is the identity");

    state.loss_channel(0, 0.1, &mut rng());
    state.x(0);

    let overlap = state.trace(&zero_state());
    assert!(
        (overlap + 0.8).abs() < 1e-10,
        "expected −0.8, got {overlap}"
    );
}

/// Old's `test_ghz_final_loss`: a GHZ backward propagation with loss applied at
/// the end.
#[test]
fn ghz_final_loss() {
    let p_l = 0.1_f64;
    let mut state = sum(2, &[("ZZ", 1.0)]);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    // Identity Cliffords must not disturb loss.
    state.x(0);
    state.x(1);
    state.x(0);
    state.x(1);

    state.loss_channel(0, p_l, &mut rng());
    state.loss_channel(1, p_l, &mut rng());

    state.cnot(0, 1);
    state.h(0);

    let overlap = state.trace(&zero_state());
    let prob = 0.5 + 0.5 * ((1.0 - p_l) * (1.0 - p_l) - 2.0 * p_l * (1.0 - p_l) + p_l * p_l);
    assert!((overlap - prob).abs() < 1e-10, "got {overlap}, want {prob}");
}

/// Old's `test_ghz`: loss both before and after the entangling gate.
#[test]
fn ghz_loss_mid_circuit() {
    let p_l = 0.1_f64;
    let mut state = sum(2, &[("ZZ", 1.0)]);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);
    state.loss_channel(0, p_l, &mut rng());
    state.loss_channel(1, p_l, &mut rng());
    state.cnot(0, 1);
    state.loss_channel(0, 2.0 * p_l, &mut rng());
    state.h(0);

    let overlap = state.trace(&zero_state());
    let prob = 2.0 * p_l
        + (1.0 - 2.0 * p_l)
            * (0.5 + 0.5 * ((1.0 - p_l) * (1.0 - p_l) - 2.0 * p_l * (1.0 - p_l) + p_l * p_l));
    assert!((overlap - prob).abs() < 1e-10, "got {overlap}, want {prob}");
}

// --- MaxLossWeight (old `test_loss_truncation`) -----------------------------

#[test]
fn max_loss_weight_truncation() {
    let mut state: LossyPauliSum<f64, MaxLossWeight> = Sum::with_policy(3, MaxLossWeight(2));
    state += (lw("ZZZ"), 1.0);
    for q in 0..3 {
        state.reset_loss_channel(q);
    }
    for q in 0..3 {
        state.loss_channel(q, 0.1, &mut rng());
    }

    let original_len = state.len();
    state.truncate();
    assert_eq!(
        state.len(),
        original_len - 1,
        "exactly the loss-weight-3 term is dropped"
    );
}

#[test]
fn max_loss_weight_defaults_and_sentinel() {
    // Old's default is 10, not `usize::MAX`.
    assert_eq!(MaxLossWeight::default(), MaxLossWeight(10));
    assert_eq!(MaxLossWeight::default().max_loss_weight(), 10);

    // `usize::MAX` disables the pass — nothing is dropped even at high loss
    // weight.
    let mut state: LossyPauliSum<f64, MaxLossWeight> =
        Sum::with_policy(3, MaxLossWeight(usize::MAX));
    state += (lw("LLL"), 1.0);
    state.truncate();
    assert_eq!(state.len(), 1);
}

// --- correlated_loss_channel ------------------------------------------------

/// The both-present arm is a pure in-place scale (`1 − 2p₁ − p₀`) that can
/// neither insert nor remove.
#[test]
fn correlated_loss_both_present_is_a_pure_scale() {
    let mut state = sum(2, &[("ZZ", 1.0), ("XI", 2.0)]);
    state.correlated_loss_channel(0, 1, [0.01, 0.02, 0.03], &mut rng());
    let f = 1.0 - 2.0 * 0.02 - 0.01;
    assert_eq!(state.len(), 2);
    assert!((state.get(&lw("ZZ")).unwrap() - f).abs() < 1e-12);
    assert!((state.get(&lw("XI")).unwrap() - 2.0 * f).abs() < 1e-12);
}

/// The one-already-lost arm: a single branch weighted `p[1]`, survivor scaled by
/// `1 − p[2]` (old's pairing, flagged as suspected old bug 4 and reproduced).
#[test]
fn correlated_loss_one_already_lost() {
    let mut state = sum(2, &[("LZ", 1.0)]);
    state.correlated_loss_channel(0, 1, [0.01, 0.02, 0.03], &mut rng());
    assert_eq!(state.len(), 2);
    assert!((state.get(&lw("LZ")).unwrap() - (1.0 - 0.03)).abs() < 1e-12);
    assert!((state.get(&lw("IZ")).unwrap() - 0.02).abs() < 1e-12);
}

/// The both-lost arm emits **three** branches — the crate's only multi-branch
/// channel, and the sole driver of the size-directed merge — and leaves the
/// survivor unscaled (old).
#[test]
fn correlated_loss_both_lost_emits_three_branches() {
    let mut state = sum(2, &[("LL", 1.0)]);
    state.correlated_loss_channel(0, 1, [0.5, 0.02, 0.25], &mut rng());
    assert_eq!(state.len(), 4, "LL + IL + LI + II");
    assert!((state.get(&lw("LL")).unwrap() - 1.0).abs() < 1e-12);
    assert!((state.get(&lw("IL")).unwrap() - 0.25).abs() < 1e-12);
    assert!((state.get(&lw("LI")).unwrap() - 0.25).abs() < 1e-12);
    assert!((state.get(&lw("II")).unwrap() - 0.5).abs() < 1e-12);
}

/// The multi-branch merge must accumulate onto an existing key, not replace it.
#[test]
fn correlated_loss_branches_accumulate() {
    let mut state = sum(2, &[("LL", 1.0), ("II", 3.0)]);
    state.correlated_loss_channel(0, 1, [0.5, 0.0, 0.25], &mut rng());
    // `II` gets 3.0 (unchanged by the both-present arm at p = [0.5, 0, 0.25]?)
    // — no: `II` is both-present, so it is scaled by 1 − 0 − 0.5 = 0.5 → 1.5,
    // and then the `LL` term's double-recovery branch adds 0.5.
    assert!((state.get(&lw("II")).unwrap() - (1.5 + 0.5)).abs() < 1e-12);
}

// --- Feature 11: the same kernels serve the lossy word ----------------------

/// The rotation kernels run on the lossy key too, skipping lost sites (old's
/// `get_lbit` early-out).
#[test]
fn rotations_skip_lost_sites() {
    let mut state = sum(2, &[("LZ", 1.0)]);
    state.rx(0, 0.7); // qubit 0 is lost → no-op
    assert_eq!(state.len(), 1);
    assert!((state.get(&lw("LZ")).unwrap() - 1.0).abs() < 1e-12);

    state.rx(1, 0.7); // qubit 1 carries Z → branches to Y
    assert_eq!(state.len(), 2);
    assert!((state.get(&lw("LZ")).unwrap() - 0.7_f64.cos()).abs() < 1e-12);
    assert!((state.get(&lw("LY")).unwrap() - 0.7_f64.sin()).abs() < 1e-12);
}
