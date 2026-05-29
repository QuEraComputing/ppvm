//! Statistical equivalence tests between the noise-aware `Sampler` (sum
//! backend) and the pure `GeneralizedTableau` simulator.
//!
//! The sum backend evolves a probability-weighted collection of branches
//! through `loss_channel`, `pauli_error`, and `depolarize` (one branch
//! per error outcome) and then samples a branch + measurement pair per shot.
//! The pure tableau applies each channel stochastically inside a single trajectory.
//! In the limit of many shots both must yield the same joint distribution
//! over `Vec<Option<bool>>` measurement outcomes.
//!
//! Each test runs N shots on both backends with deterministic seeds and
//! checks that the total variation distance between the empirical
//! distributions is below a finite-sample threshold.

use std::collections::{HashMap, HashSet};

use ppvm_runtime::config::fxhash::ByteF64;
use ppvm_runtime::traits::{
    Clifford, CliffordExtensions, CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel,
    LossyMeasure, PauliError, Reset, RotationOne, RotationTwo, TGate, TwoQubitPauliError, U3Gate,
};
use ppvm_tableau::measure_all::LossyMeasureAll;
use ppvm_tableau::prelude::*;
use ppvm_tableau_sum::prelude::*;

type Cfg = ByteF64<1>;
type TabSum = GeneralizedTableauSum<Cfg, u128>;
type Tab = GeneralizedTableau<Cfg, u128>;

const SEED_SUM: u64 = 0xC0FFEE_u64;
const SEED_PURE: u64 = 0xDEADBEEF_u64;

fn frequencies(shots: &[Vec<Option<bool>>]) -> HashMap<Vec<Option<bool>>, f64> {
    let n = shots.len() as f64;
    let mut m: HashMap<Vec<Option<bool>>, f64> = HashMap::new();
    for s in shots {
        *m.entry(s.clone()).or_insert(0.0) += 1.0 / n;
    }
    m
}

fn tvd(a: &HashMap<Vec<Option<bool>>, f64>, b: &HashMap<Vec<Option<bool>>, f64>) -> f64 {
    let mut keys: HashSet<&Vec<Option<bool>>> = HashSet::new();
    keys.extend(a.keys());
    keys.extend(b.keys());
    let mut total = 0.0;
    for k in keys {
        let av = a.get(k).copied().unwrap_or(0.0);
        let bv = b.get(k).copied().unwrap_or(0.0);
        total += (av - bv).abs();
    }
    0.5 * total
}

fn run_sum<F>(
    n_qubits: usize,
    shots: usize,
    sum_cutoff: f64,
    mut circuit: F,
) -> Vec<Vec<Option<bool>>>
where
    F: FnMut(&mut TabSum),
{
    let mut tab: TabSum =
        GeneralizedTableauSum::new_with_seed(n_qubits, 1e-12, sum_cutoff, SEED_SUM);
    circuit(&mut tab);
    tab.sampler().sample_shots(shots)
}

fn run_pure<F>(n_qubits: usize, shots: usize, mut circuit: F) -> Vec<Vec<Option<bool>>>
where
    F: FnMut(&mut Tab),
{
    (0..shots as u64)
        .map(|i| {
            let mut t: Tab =
                GeneralizedTableau::new_with_seed(n_qubits, 1e-12, SEED_PURE.wrapping_add(i));
            circuit(&mut t);
            t.measure_all()
        })
        .collect()
}

/// Assert TVD < `tol`. On failure dump both distributions for debugging.
#[track_caller]
fn assert_distributions_match(
    sum: &[Vec<Option<bool>>],
    pure: &[Vec<Option<bool>>],
    tol: f64,
    label: &str,
) {
    let fs = frequencies(sum);
    let fp = frequencies(pure);
    let d = tvd(&fs, &fp);
    if d >= tol {
        let mut keys: Vec<_> = fs
            .keys()
            .chain(fp.keys())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        keys.sort();
        let mut report = String::new();
        for k in keys {
            report.push_str(&format!(
                "  {:?}  sum={:.4}  pure={:.4}\n",
                k,
                fs.get(k).copied().unwrap_or(0.0),
                fp.get(k).copied().unwrap_or(0.0)
            ));
        }
        panic!("[{label}] TVD = {d:.4} >= tol {tol}\n{report}");
    }
}

// ---------------------------------------------------------------------------
// Single-qubit loss
// ---------------------------------------------------------------------------

#[test]
fn loss_channel_after_hadamard() {
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.loss_channel(0, 0.3);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.loss_channel(0, 0.3);
    });
    // Three outcomes: None, Some(false), Some(true). 5σ ≈ 0.025 per bin.
    assert_distributions_match(&sum, &pure, 0.04, "loss_channel_after_hadamard");
}

#[test]
fn loss_channel_p_one_marks_all_lost() {
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.x(0);
        t.loss_channel(0, 1.0);
    });
    let pure = run_pure(1, shots, |t| {
        t.x(0);
        t.loss_channel(0, 1.0);
    });
    // Deterministic: every shot must be None.
    assert!(sum.iter().all(|s| s[0].is_none()));
    assert!(pure.iter().all(|s| s[0].is_none()));
}

#[test]
fn loss_channel_p_zero_no_loss() {
    let shots = 2000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.loss_channel(0, 0.0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.loss_channel(0, 0.0);
    });
    assert!(sum.iter().all(|s| s[0].is_some()));
    assert!(pure.iter().all(|s| s[0].is_some()));
    assert_distributions_match(&sum, &pure, 0.05, "loss_channel_p_zero_no_loss");
}

// ---------------------------------------------------------------------------
// Bell pair + loss
// ---------------------------------------------------------------------------

#[test]
fn bell_pair_with_loss_on_q0() {
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.loss_channel(0, 0.3);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.loss_channel(0, 0.3);
    });
    // 4 possible outcomes: (None,0), (None,1), (0,0), (1,1). 5σ ≈ 0.03 per bin.
    assert_distributions_match(&sum, &pure, 0.05, "bell_pair_with_loss_on_q0");
}

#[test]
fn bell_pair_with_loss_on_both_qubits() {
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.loss_channel(0, 0.2);
        t.loss_channel(1, 0.2);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.loss_channel(0, 0.2);
        t.loss_channel(1, 0.2);
    });
    assert_distributions_match(&sum, &pure, 0.06, "bell_pair_with_loss_on_both_qubits");
}

// ---------------------------------------------------------------------------
// Single-qubit depolarizing channel
// ---------------------------------------------------------------------------

#[test]
fn depolarize_on_ground_state() {
    let shots = 8000;
    let p = 0.6_f64;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.depolarize(0, p);
    });
    let pure = run_pure(1, shots, |t| {
        t.depolarize(0, p);
    });
    // P(1) should converge to 2p/3 on both sides.
    let ones_sum = sum.iter().filter(|s| s[0] == Some(true)).count() as f64 / shots as f64;
    let ones_pure = pure.iter().filter(|s| s[0] == Some(true)).count() as f64 / shots as f64;
    let expected = 2.0 * p / 3.0;
    // 5σ for p≈0.4, N=8000 ≈ 0.027
    assert!(
        (ones_sum - expected).abs() < 0.04,
        "sum P(1)={ones_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (ones_pure - expected).abs() < 0.04,
        "pure P(1)={ones_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(&sum, &pure, 0.04, "depolarize_on_ground_state");
}

#[test]
fn depolarize_on_plus_state() {
    let shots = 4000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.depolarize(0, 0.5);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.depolarize(0, 0.5);
    });
    // |+⟩ is invariant under depolarization in the Z basis: P(0)=P(1)=0.5.
    assert_distributions_match(&sum, &pure, 0.05, "depolarize_on_plus_state");
}

// ---------------------------------------------------------------------------
// Single-qubit Pauli error channel
// ---------------------------------------------------------------------------

#[test]
fn pauli_error_on_ground_state_nonuniform() {
    let shots = 8000;
    let p = [0.15_f64, 0.25_f64, 0.35_f64];
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.pauli_error(0, p);
    });
    let pure = run_pure(1, shots, |t| {
        t.pauli_error(0, p);
    });

    let ones_sum = sum.iter().filter(|s| s[0] == Some(true)).count() as f64 / shots as f64;
    let ones_pure = pure.iter().filter(|s| s[0] == Some(true)).count() as f64 / shots as f64;
    let expected = p[0] + p[1];
    assert!(
        (ones_sum - expected).abs() < 0.04,
        "sum P(1)={ones_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (ones_pure - expected).abs() < 0.04,
        "pure P(1)={ones_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(&sum, &pure, 0.04, "pauli_error_on_ground_state_nonuniform");
}

#[test]
fn pauli_error_on_lost_qubit_is_noop() {
    let shots = 1000;
    let p = [0.2_f64, 0.3_f64, 0.1_f64];
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.pauli_error(0, p);
    });
    let pure = run_pure(1, shots, |t| {
        t.loss_channel(0, 1.0);
        t.pauli_error(0, p);
    });

    assert!(sum.iter().all(|s| s[0].is_none()));
    assert!(pure.iter().all(|s| s[0].is_none()));
}

// ---------------------------------------------------------------------------
// Bell pair + depolarizing channel
// ---------------------------------------------------------------------------

#[test]
fn bell_pair_with_depolarize_on_q0() {
    // After depolarize(q0, p) on a Bell pair (|00⟩+|11⟩)/√2:
    //   I, Z keep correlation (|00⟩+|11⟩ or |00⟩-|11⟩) → measurements agree
    //   X, Y break correlation → measurements disagree
    // So P(same) = 1 - 2p/3 and P(diff) = 2p/3.
    let shots = 8000;
    let p = 0.3_f64;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.depolarize(0, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.depolarize(0, p);
    });
    let same_sum = sum.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let same_pure = pure.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let expected = 1.0 - 2.0 * p / 3.0;
    assert!(
        (same_sum - expected).abs() < 0.04,
        "sum P(same)={same_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (same_pure - expected).abs() < 0.04,
        "pure P(same)={same_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(&sum, &pure, 0.05, "bell_pair_with_depolarize_on_q0");
}

#[test]
fn bell_pair_with_pauli_error_on_q0_nonuniform() {
    let shots = 8000;
    let p = [0.1_f64, 0.2_f64, 0.3_f64];
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.pauli_error(0, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.pauli_error(0, p);
    });

    let same_sum = sum.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let same_pure = pure.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let expected = 1.0 - p[0] - p[1];
    assert!(
        (same_sum - expected).abs() < 0.04,
        "sum P(same)={same_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (same_pure - expected).abs() < 0.04,
        "pure P(same)={same_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.05,
        "bell_pair_with_pauli_error_on_q0_nonuniform",
    );
}

// ---------------------------------------------------------------------------
// Mixed loss + depolarize
// ---------------------------------------------------------------------------

#[test]
fn loss_then_depolarize_three_qubits() {
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        t.loss_channel(0, 0.2);
        t.depolarize(1, 0.15);
        t.depolarize(2, 0.15);
    });
    let pure = run_pure(3, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        t.loss_channel(0, 0.2);
        t.depolarize(1, 0.15);
        t.depolarize(2, 0.15);
    });
    // Outcome space: 3 values per qubit ⇒ up to 27 bins; many will have
    // ~0 mass so per-bin error is dominated by the high-mass ones.
    assert_distributions_match(&sum, &pure, 0.08, "loss_then_depolarize_three_qubits");
}

#[test]
fn ghz_three_qubits_with_per_qubit_noise() {
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(0, 2);
        for q in 0..3 {
            t.depolarize(q, 0.1);
            t.loss_channel(q, 0.05);
        }
    });
    let pure = run_pure(3, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(0, 2);
        for q in 0..3 {
            t.depolarize(q, 0.1);
            t.loss_channel(q, 0.05);
        }
    });
    assert_distributions_match(&sum, &pure, 0.08, "ghz_three_qubits_with_per_qubit_noise");
}

#[test]
fn repeated_depolarize_creates_many_branches() {
    // Apply depolarize repeatedly on the same qubit. The sum's branch
    // deduplication must keep statistics correct under accumulated noise.
    let shots = 8000;
    let p = 0.1_f64;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        for _ in 0..5 {
            t.depolarize(0, p);
            t.depolarize(1, p);
        }
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        for _ in 0..5 {
            t.depolarize(0, p);
            t.depolarize(1, p);
        }
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.06,
        "repeated_depolarize_creates_many_branches",
    );
}

#[test]
fn clifford_layer_with_sqrt_gates_and_noise() {
    // Exercise the CliffordExtensions (sqrt_x, sqrt_y, ...) wired through
    // the sum backend together with both noise channels.
    let shots = 8000;
    let sum = run_sum(4, shots, 1e-12, |t| {
        t.sqrt_y(0);
        t.sqrt_y(1);
        t.sqrt_y(2);
        t.sqrt_y(3);
        t.cz(0, 1);
        t.cz(2, 3);
        for q in 0..4 {
            t.depolarize(q, 0.08);
        }
        t.cz(1, 2);
        t.sqrt_x_adj(0);
        t.sqrt_x_adj(3);
        for q in 0..4 {
            t.loss_channel(q, 0.05);
        }
    });
    let pure = run_pure(4, shots, |t| {
        t.sqrt_y(0);
        t.sqrt_y(1);
        t.sqrt_y(2);
        t.sqrt_y(3);
        t.cz(0, 1);
        t.cz(2, 3);
        for q in 0..4 {
            t.depolarize(q, 0.08);
        }
        t.cz(1, 2);
        t.sqrt_x_adj(0);
        t.sqrt_x_adj(3);
        for q in 0..4 {
            t.loss_channel(q, 0.05);
        }
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.08,
        "clifford_layer_with_sqrt_gates_and_noise",
    );
}

// ---------------------------------------------------------------------------
// Reset
// ---------------------------------------------------------------------------

#[test]
#[ignore = "TODO"]
fn reset_after_hadamard_collapses_to_zero() {
    // Reset is `measure + flip-if-1`, so after any single-qubit state it
    // forces the qubit to |0⟩. Both backends must give Some(false) every shot.
    let shots = 2000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.reset(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.reset(0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(false)));
    assert!(pure.iter().all(|s| s[0] == Some(false)));
}

// Known failure: `GeneralizedTableauSum::reset` delegates to per-entry
// `GeneralizedTableau::reset`, which is `measure + flip-if-1`. The inner
// measurement consumes the entry's RNG to commit to one outcome, so for any
// state where q0 is entangled with other qubits the "build once, sample many"
// model collapses the entry to a single product state and the other branch's
// probability mass is lost.
//
// Pure backend: rebuilds the state per shot, so the stochastic reset is
// spread across shots — 50/50 between (0,0) and (0,1).
// Sum backend:  one entry; reset picks one outcome; all 8000 shots return it.
//
// Correct fix is to make `reset` branch on the sum (each entry → 2 sum-level
// entries with probabilities ⟨0|ρ|0⟩ and ⟨1|ρ|1⟩, both with q0=0 post-state),
// analogous to how `loss_channel` / `depolarize` already work. Tracking this
// here so future work on `Reset` lands together with re-enabling the test.
#[test]
#[ignore = "Reset on the sum backend doesn't branch; see comment above"]
fn reset_bell_pair_q0_decorrelates() {
    // Bell pair then reset(q0): q0 is forced to 0 and the reset's
    // measurement collapses q1 to a definite (but random) value. So
    // q0 = 0 always, q1 is 50/50, and the two are independent in the
    // ensemble.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.reset(0);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.reset(0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(false)));
    assert!(pure.iter().all(|s| s[0] == Some(false)));
    assert_distributions_match(&sum, &pure, 0.04, "reset_bell_pair_q0_decorrelates");
}

#[test]
#[ignore = "TODO"]
fn reset_after_depolarize_still_zero() {
    // Depolarize then reset: regardless of which Pauli error fired,
    // reset projects back to |0⟩ deterministically.
    let shots = 4000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.depolarize(0, 0.5);
        t.reset(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.depolarize(0, 0.5);
        t.reset(0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(false)));
    assert!(pure.iter().all(|s| s[0] == Some(false)));
}

#[test]
#[ignore = "TODO"]
fn reset_on_lost_qubit_is_no_op() {
    // GeneralizedTableau::reset is `measure + flip`; measurement on a lost
    // qubit returns None, so the reset leaves is_lost set. Both backends
    // must still report None for that qubit.
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.reset(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.loss_channel(0, 1.0);
        t.reset(0);
    });
    assert!(sum.iter().all(|s| s[0].is_none()));
    assert!(pure.iter().all(|s| s[0].is_none()));
}

// ---------------------------------------------------------------------------
// Single-qubit rotations (RotationOne)
// ---------------------------------------------------------------------------

#[test]
fn rx_pi_flips_zero_to_one() {
    // RX(π) on |0⟩ → |1⟩ deterministically. No branching inside the gate.
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.rx(0, std::f64::consts::PI);
    });
    let pure = run_pure(1, shots, |t| {
        t.rx(0, std::f64::consts::PI);
    });
    assert!(sum.iter().all(|s| s[0] == Some(true)));
    assert!(pure.iter().all(|s| s[0] == Some(true)));
}

#[test]
fn rx_half_pi_creates_unbiased_distribution() {
    // RX(π/2)|0⟩ has equal |0⟩/|1⟩ Z-basis probabilities.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.rx(0, std::f64::consts::FRAC_PI_2);
    });
    let pure = run_pure(1, shots, |t| {
        t.rx(0, std::f64::consts::FRAC_PI_2);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "rx_half_pi_creates_unbiased_distribution",
    );
}

#[test]
fn rz_alone_does_not_change_z_basis_distribution() {
    // RZ is diagonal in the Z basis: |0⟩ stays at P(0)=1.
    let shots = 2000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.rz(0, 0.7);
    });
    let pure = run_pure(1, shots, |t| {
        t.rz(0, 0.7);
    });
    assert!(sum.iter().all(|s| s[0] == Some(false)));
    assert!(pure.iter().all(|s| s[0] == Some(false)));
}

#[test]
fn ry_then_rz_then_ry_with_depolarize() {
    // Non-trivial single-qubit angle sequence + depolarize to check that
    // RotationOne composes with noise the same way on both backends.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.ry(0, 0.41 * std::f64::consts::PI);
        t.rz(0, 0.23 * std::f64::consts::PI);
        t.ry(0, 0.17 * std::f64::consts::PI);
        t.depolarize(0, 0.12);
    });
    let pure = run_pure(1, shots, |t| {
        t.ry(0, 0.41 * std::f64::consts::PI);
        t.rz(0, 0.23 * std::f64::consts::PI);
        t.ry(0, 0.17 * std::f64::consts::PI);
        t.depolarize(0, 0.12);
    });
    assert_distributions_match(&sum, &pure, 0.04, "ry_then_rz_then_ry_with_depolarize");
}

#[test]
fn rx_then_loss_two_qubits() {
    // Single-qubit rotation followed by loss on one of two qubits; verifies
    // RotationOne plays well with the sum's loss branches.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.rx(0, std::f64::consts::FRAC_PI_2);
        t.cnot(0, 1);
        t.loss_channel(0, 0.25);
    });
    let pure = run_pure(2, shots, |t| {
        t.rx(0, std::f64::consts::FRAC_PI_2);
        t.cnot(0, 1);
        t.loss_channel(0, 0.25);
    });
    assert_distributions_match(&sum, &pure, 0.05, "rx_then_loss_two_qubits");
}

// ---------------------------------------------------------------------------
// Two-qubit rotations (RotationTwo)
// ---------------------------------------------------------------------------

#[test]
fn rxx_pi_flips_both_qubits() {
    // RXX(π)|00⟩ = -i|11⟩: deterministic flip.
    let shots = 1000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.rxx(0, 1, std::f64::consts::PI);
    });
    let pure = run_pure(2, shots, |t| {
        t.rxx(0, 1, std::f64::consts::PI);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(true), Some(true)]));
    assert!(pure.iter().all(|s| s == &vec![Some(true), Some(true)]));
}

#[test]
fn rxx_half_pi_correlated_outcomes() {
    // RXX(π/2)|00⟩ = (|00⟩ - i|11⟩)/√2: q0 and q1 always agree.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.rxx(0, 1, std::f64::consts::FRAC_PI_2);
    });
    let pure = run_pure(2, shots, |t| {
        t.rxx(0, 1, std::f64::consts::FRAC_PI_2);
    });
    assert!(sum.iter().all(|s| s[0] == s[1]));
    assert!(pure.iter().all(|s| s[0] == s[1]));
    assert_distributions_match(&sum, &pure, 0.04, "rxx_half_pi_correlated_outcomes");
}

#[test]
fn ryy_half_pi_correlated_outcomes() {
    // RYY(π/2)|00⟩ = (|00⟩ + i|11⟩)/√2: q0 and q1 always agree.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.ryy(0, 1, std::f64::consts::FRAC_PI_2);
    });
    let pure = run_pure(2, shots, |t| {
        t.ryy(0, 1, std::f64::consts::FRAC_PI_2);
    });
    assert!(sum.iter().all(|s| s[0] == s[1]));
    assert!(pure.iter().all(|s| s[0] == s[1]));
    assert_distributions_match(&sum, &pure, 0.04, "ryy_half_pi_correlated_outcomes");
}

#[test]
fn rzz_diagonal_on_comp_basis() {
    // RZZ is diagonal in the Z basis: |00⟩ stays a pure Z eigenstate.
    let shots = 2000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.rzz(0, 1, 0.4 * std::f64::consts::PI);
    });
    let pure = run_pure(2, shots, |t| {
        t.rzz(0, 1, 0.4 * std::f64::consts::PI);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(false), Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![Some(false), Some(false)]));
}

#[test]
fn rxx_with_depolarize_breaks_correlation() {
    // RXX(π/2)|00⟩ correlates q0 and q1. Depolarize(q0, p) restores the
    // 2p/3 mismatch probability. Compare sum vs pure across all 4 bins.
    let shots = 8000;
    let p = 0.2_f64;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.rxx(0, 1, std::f64::consts::FRAC_PI_2);
        t.depolarize(0, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.rxx(0, 1, std::f64::consts::FRAC_PI_2);
        t.depolarize(0, p);
    });
    let agree_sum = sum.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let agree_pure = pure.iter().filter(|s| s[0] == s[1]).count() as f64 / shots as f64;
    let expected = 1.0 - 2.0 * p / 3.0;
    assert!(
        (agree_sum - expected).abs() < 0.04,
        "sum P(agree)={agree_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (agree_pure - expected).abs() < 0.04,
        "pure P(agree)={agree_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(&sum, &pure, 0.05, "rxx_with_depolarize_breaks_correlation");
}

#[test]
fn rxy_three_qubit_chain_with_loss() {
    // Mixes the rxy/ryz axes (less commonly exercised) with a loss channel
    // to ensure the generic rotate_2 path is consistent.
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.rxy(0, 1, 0.3 * std::f64::consts::PI);
        t.ryz(1, 2, 0.25 * std::f64::consts::PI);
        t.loss_channel(2, 0.1);
    });
    let pure = run_pure(3, shots, |t| {
        t.rxy(0, 1, 0.3 * std::f64::consts::PI);
        t.ryz(1, 2, 0.25 * std::f64::consts::PI);
        t.loss_channel(2, 0.1);
    });
    assert_distributions_match(&sum, &pure, 0.08, "rxy_three_qubit_chain_with_loss");
}

// ---------------------------------------------------------------------------
// Two-qubit Pauli error channel
// ---------------------------------------------------------------------------
//
// Probability-array index layout (matches GeneralizedTableau and the sum
// backend): IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ.

#[test]
fn two_qubit_pauli_error_zero_prob_is_noop() {
    let shots = 1000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, [0.0; 15]);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, [0.0; 15]);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(false), Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![Some(false), Some(false)]));
}

#[test]
fn two_qubit_pauli_error_xx_certain_flips_both() {
    // p[4] = XX with probability 1.0: |00⟩ → |11⟩ deterministically.
    let shots = 1000;
    let mut p = [0.0_f64; 15];
    p[4] = 1.0;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(true), Some(true)]));
    assert!(pure.iter().all(|s| s == &vec![Some(true), Some(true)]));
}

#[test]
fn two_qubit_pauli_error_zz_invariant_in_z_basis() {
    // p[14] = ZZ with probability 1.0: |00⟩ is a ZZ eigenstate, no flips.
    let shots = 1000;
    let mut p = [0.0_f64; 15];
    p[14] = 1.0;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(false), Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![Some(false), Some(false)]));
}

#[test]
fn two_qubit_pauli_error_xi_flips_q0_only() {
    // p[3] = XI with probability 1.0: q0 flips, q1 stays.
    let shots = 1000;
    let mut p = [0.0_f64; 15];
    p[3] = 1.0;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert!(sum.iter().all(|s| s == &vec![Some(true), Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![Some(true), Some(false)]));
}

#[test]
fn two_qubit_pauli_error_uniform_on_ground_state() {
    // Uniform 1/15 weighting: total probability 1, every non-identity Pauli
    // equally likely. Flip distribution on |00⟩:
    //   |00⟩: I⊗Z, Z⊗I, Z⊗Z       (3/15)
    //   |01⟩: I⊗{X,Y}, Z⊗{X,Y}    (4/15)
    //   |10⟩: {X,Y}⊗I, {X,Y}⊗Z    (4/15)
    //   |11⟩: {X,Y}⊗{X,Y}         (4/15)
    let shots = 8000;
    let p = [1.0 / 15.0_f64; 15];
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "two_qubit_pauli_error_uniform_on_ground_state",
    );
}

#[test]
fn two_qubit_pauli_error_nonuniform_on_ground_state() {
    // Nontrivial heterogeneous probabilities: exercises that each p[k] lands
    // at the right Pauli pair on both backends.
    let shots = 8000;
    let p = [
        0.01, 0.02, 0.03, // IX, IY, IZ
        0.04, 0.05, 0.06, 0.07, // XI, XX, XY, XZ
        0.05, 0.04, 0.03, 0.02, // YI, YX, YY, YZ
        0.01, 0.02, 0.03, 0.04, // ZI, ZX, ZY, ZZ
    ];
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "two_qubit_pauli_error_nonuniform_on_ground_state",
    );
}

#[test]
fn two_qubit_pauli_error_on_lost_qubit_is_noop() {
    // If either input qubit is lost the channel short-circuits to a no-op.
    let shots = 1000;
    let p = [1.0 / 15.0_f64; 15];
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.loss_channel(0, 1.0);
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert!(sum.iter().all(|s| s == &vec![None, Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![None, Some(false)]));
}

#[test]
fn depolarize2_on_ground_state() {
    // depolarize2(p) = two_qubit_pauli_error with uniform p/15 weights. From
    // |00⟩, 12 of the 15 non-identity Paulis cause a Z-basis flip on at
    // least one qubit (only IZ, ZI, ZZ leave |00⟩ invariant), so
    // P(any flip) = 12p/15 = 4p/5.
    let shots = 8000;
    let p = 0.6_f64;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.depolarize2(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.depolarize2(0, 1, p);
    });
    let flipped_sum =
        sum.iter().filter(|s| s != &&vec![Some(false), Some(false)]).count() as f64 / shots as f64;
    let flipped_pure = pure
        .iter()
        .filter(|s| s != &&vec![Some(false), Some(false)])
        .count() as f64
        / shots as f64;
    let expected = 4.0 * p / 5.0;
    assert!(
        (flipped_sum - expected).abs() < 0.04,
        "sum P(flip)={flipped_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (flipped_pure - expected).abs() < 0.04,
        "pure P(flip)={flipped_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(&sum, &pure, 0.04, "depolarize2_on_ground_state");
}

#[test]
fn bell_pair_with_two_qubit_pauli_error_nonuniform() {
    // Bell pair + nontrivial two-qubit Pauli noise; statistics must agree.
    let shots = 8000;
    let p = [
        0.02, 0.03, 0.04, // IX, IY, IZ
        0.05, 0.06, 0.04, 0.03, // XI, XX, XY, XZ
        0.04, 0.05, 0.03, 0.02, // YI, YX, YY, YZ
        0.02, 0.03, 0.04, 0.05, // ZI, ZX, ZY, ZZ
    ];
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.two_qubit_pauli_error(0, 1, p);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.two_qubit_pauli_error(0, 1, p);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.05,
        "bell_pair_with_two_qubit_pauli_error_nonuniform",
    );
}

// ---------------------------------------------------------------------------
// Correlated two-qubit loss channel
// ---------------------------------------------------------------------------
//
// Probability layout (matches the trait spec and GeneralizedTableau):
//   p[0] = P(lose both | both present)
//   p[1] = P(lose either one | both present)      → split 50/50 across q0/q1
//   p[2] = P(lose remaining | the other was lost prior)

#[test]
fn correlated_loss_channel_zero_prob_is_noop() {
    let shots = 1000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0]);
    });
    let pure = run_pure(2, shots, |t| {
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0]);
    });
    // Both qubits remain in the qubit subspace and stay |0⟩.
    assert!(sum.iter().all(|s| s == &vec![Some(false), Some(false)]));
    assert!(pure.iter().all(|s| s == &vec![Some(false), Some(false)]));
}

#[test]
fn correlated_loss_channel_both_certain() {
    // p[0] = 1.0: both qubits deterministically lost.
    let shots = 1000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0]);
    });
    let pure = run_pure(2, shots, |t| {
        t.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0]);
    });
    assert!(sum.iter().all(|s| s == &vec![None, None]));
    assert!(pure.iter().all(|s| s == &vec![None, None]));
}

#[test]
fn correlated_loss_channel_single_loss_certain_is_5050_between_qubits() {
    // p[1] = 1.0: exactly one of the two qubits is lost, 50/50 which one.
    // Never both, never neither.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]);
    });
    let pure = run_pure(2, shots, |t| {
        t.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]);
    });
    // Sanity: no shot has both lost or neither lost.
    assert!(
        sum.iter().all(|s| s[0].is_none() ^ s[1].is_none()),
        "sum: every shot must lose exactly one qubit"
    );
    assert!(
        pure.iter().all(|s| s[0].is_none() ^ s[1].is_none()),
        "pure: every shot must lose exactly one qubit"
    );
    // q0-lost frequency should converge to 1/2 on both backends.
    let q0_lost_sum = sum.iter().filter(|s| s[0].is_none()).count() as f64 / shots as f64;
    let q0_lost_pure = pure.iter().filter(|s| s[0].is_none()).count() as f64 / shots as f64;
    assert!(
        (q0_lost_sum - 0.5).abs() < 0.04,
        "sum P(q0 lost)={q0_lost_sum:.4}"
    );
    assert!(
        (q0_lost_pure - 0.5).abs() < 0.04,
        "pure P(q0 lost)={q0_lost_pure:.4}"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "correlated_loss_channel_single_loss_certain_is_5050_between_qubits",
    );
}

#[test]
fn correlated_loss_channel_marginals_on_ground_state() {
    // Mid-probability values, |00⟩ input. Check both the analytic outcome
    // probabilities and TVD against pure.
    // p[0] = 0.20, p[1] = 0.40, p[2] is unused (no qubit pre-lost).
    // Expected outcome probabilities:
    //   P(both lost)    = 0.20
    //   P(only q0 lost) = 0.20  (p[1]/2)
    //   P(only q1 lost) = 0.20
    //   P(none lost)    = 0.40
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.correlated_loss_channel(0, 1, [0.20, 0.40, 0.0]);
    });
    let pure = run_pure(2, shots, |t| {
        t.correlated_loss_channel(0, 1, [0.20, 0.40, 0.0]);
    });
    for (label, data) in [("sum", &sum), ("pure", &pure)] {
        let both = data.iter().filter(|s| s[0].is_none() && s[1].is_none()).count() as f64
            / shots as f64;
        let only_q0 = data.iter().filter(|s| s[0].is_none() && s[1].is_some()).count() as f64
            / shots as f64;
        let only_q1 = data.iter().filter(|s| s[0].is_some() && s[1].is_none()).count() as f64
            / shots as f64;
        let none = data.iter().filter(|s| s[0].is_some() && s[1].is_some()).count() as f64
            / shots as f64;
        assert!((both - 0.20).abs() < 0.04, "{label} P(both)={both:.4}");
        assert!(
            (only_q0 - 0.20).abs() < 0.04,
            "{label} P(only q0)={only_q0:.4}"
        );
        assert!(
            (only_q1 - 0.20).abs() < 0.04,
            "{label} P(only q1)={only_q1:.4}"
        );
        assert!((none - 0.40).abs() < 0.04, "{label} P(none)={none:.4}");
    }
    assert_distributions_match(
        &sum,
        &pure,
        0.05,
        "correlated_loss_channel_marginals_on_ground_state",
    );
}

#[test]
fn correlated_loss_channel_preexisting_loss_falls_back_to_p2() {
    // q0 is already lost (prior loss_channel with p=1.0). Channel must
    // reduce to a single-qubit Bernoulli loss on q1 with probability p[2].
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.correlated_loss_channel(0, 1, [0.5, 0.5, 0.3]);
    });
    let pure = run_pure(2, shots, |t| {
        t.loss_channel(0, 1.0);
        t.correlated_loss_channel(0, 1, [0.5, 0.5, 0.3]);
    });
    // q0 must be lost every shot. q1 must be lost with frequency p[2].
    assert!(sum.iter().all(|s| s[0].is_none()));
    assert!(pure.iter().all(|s| s[0].is_none()));
    let q1_lost_sum = sum.iter().filter(|s| s[1].is_none()).count() as f64 / shots as f64;
    let q1_lost_pure = pure.iter().filter(|s| s[1].is_none()).count() as f64 / shots as f64;
    assert!(
        (q1_lost_sum - 0.3).abs() < 0.04,
        "sum P(q1 lost)={q1_lost_sum:.4}, expected 0.3"
    );
    assert!(
        (q1_lost_pure - 0.3).abs() < 0.04,
        "pure P(q1 lost)={q1_lost_pure:.4}, expected 0.3"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "correlated_loss_channel_preexisting_loss_falls_back_to_p2",
    );
}

#[test]
fn correlated_loss_channel_both_preexisting_loss_is_noop() {
    // Both qubits pre-lost: channel must leave the state untouched.
    let shots = 1000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.loss_channel(1, 1.0);
        t.correlated_loss_channel(0, 1, [0.5, 0.5, 0.5]);
    });
    let pure = run_pure(2, shots, |t| {
        t.loss_channel(0, 1.0);
        t.loss_channel(1, 1.0);
        t.correlated_loss_channel(0, 1, [0.5, 0.5, 0.5]);
    });
    assert!(sum.iter().all(|s| s == &vec![None, None]));
    assert!(pure.iter().all(|s| s == &vec![None, None]));
}

#[test]
fn bell_pair_with_correlated_loss_channel() {
    // Bell pair (|00⟩+|11⟩)/√2, then correlated loss. Mid-probability values.
    // Statistics on both backends must agree.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.correlated_loss_channel(0, 1, [0.10, 0.20, 0.15]);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.correlated_loss_channel(0, 1, [0.10, 0.20, 0.15]);
    });
    assert_distributions_match(&sum, &pure, 0.05, "bell_pair_with_correlated_loss_channel");
}

// ---------------------------------------------------------------------------
// U3
// ---------------------------------------------------------------------------

#[test]
fn u3_x_gate() {
    // U3(π, 0, 0) acts as RY(π) on |0⟩ → |1⟩ (up to global phase).
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.u3(0, std::f64::consts::PI, 0.0, 0.0);
    });
    let pure = run_pure(1, shots, |t| {
        t.u3(0, std::f64::consts::PI, 0.0, 0.0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(true)));
    assert!(pure.iter().all(|s| s[0] == Some(true)));
}

#[test]
fn u3_random_angles_with_depolarize() {
    // Generic U3 with non-trivial (θ, φ, λ) followed by depolarize.
    let shots = 8000;
    let (theta, phi, lambda) = (
        0.37 * std::f64::consts::PI,
        0.18 * std::f64::consts::PI,
        0.51 * std::f64::consts::PI,
    );
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.u3(0, theta, phi, lambda);
        t.depolarize(0, 0.15);
    });
    let pure = run_pure(1, shots, |t| {
        t.u3(0, theta, phi, lambda);
        t.depolarize(0, 0.15);
    });
    assert_distributions_match(&sum, &pure, 0.04, "u3_random_angles_with_depolarize");
}

#[test]
fn u3_two_qubit_circuit_with_loss() {
    // U3 on q0, entangle with q1 via CNOT, loss on q1.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.u3(
            0,
            0.4 * std::f64::consts::PI,
            0.2 * std::f64::consts::PI,
            0.1 * std::f64::consts::PI,
        );
        t.cnot(0, 1);
        t.loss_channel(1, 0.2);
    });
    let pure = run_pure(2, shots, |t| {
        t.u3(
            0,
            0.4 * std::f64::consts::PI,
            0.2 * std::f64::consts::PI,
            0.1 * std::f64::consts::PI,
        );
        t.cnot(0, 1);
        t.loss_channel(1, 0.2);
    });
    assert_distributions_match(&sum, &pure, 0.06, "u3_two_qubit_circuit_with_loss");
}

// ---------------------------------------------------------------------------
// T gate
// ---------------------------------------------------------------------------

#[test]
fn t_then_h_distribution() {
    // H · T · H |0⟩ = (cos(π/8)|0⟩ - i·sin(π/8)|1⟩)·... — the resulting Z-basis
    // probability is non-trivial (≠ 0, 1/2, 1). Sum and pure must agree.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.t(0);
        t.h(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.t(0);
        t.h(0);
    });
    assert_distributions_match(&sum, &pure, 0.04, "t_then_h_distribution");
}

#[test]
fn t_h_t_adj_with_depolarize() {
    // Non-Clifford T/T† interleaved with H and depolarizing noise: verifies
    // TGate composes correctly inside the sum, end-to-end with noise.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.t(0);
        t.h(0);
        t.depolarize(0, 0.1);
        t.h(0);
        t.t_adj(0);
        t.h(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.t(0);
        t.h(0);
        t.depolarize(0, 0.1);
        t.h(0);
        t.t_adj(0);
        t.h(0);
    });
    assert_distributions_match(&sum, &pure, 0.04, "t_h_t_adj_with_depolarize");
}

// ---------------------------------------------------------------------------
// Mixed: many gate families + noise in one circuit
// ---------------------------------------------------------------------------

#[test]
#[ignore = "TODO"]
fn mixed_rotations_reset_t_noise() {
    // Exercise all newly-wired gate impls together with both noise channels.
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.rx(0, 0.3 * std::f64::consts::PI);
        t.u3(
            1,
            0.4 * std::f64::consts::PI,
            0.1 * std::f64::consts::PI,
            0.2 * std::f64::consts::PI,
        );
        t.cnot(0, 1);
        t.rxx(1, 2, 0.25 * std::f64::consts::PI);
        t.t(2);
        t.depolarize(0, 0.08);
        t.reset(0);
        t.loss_channel(2, 0.1);
    });
    let pure = run_pure(3, shots, |t| {
        t.rx(0, 0.3 * std::f64::consts::PI);
        t.u3(
            1,
            0.4 * std::f64::consts::PI,
            0.1 * std::f64::consts::PI,
            0.2 * std::f64::consts::PI,
        );
        t.cnot(0, 1);
        t.rxx(1, 2, 0.25 * std::f64::consts::PI);
        t.t(2);
        t.depolarize(0, 0.08);
        t.reset(0);
        t.loss_channel(2, 0.1);
    });
    assert_distributions_match(&sum, &pure, 0.08, "mixed_rotations_reset_t_noise");
}

#[test]
fn truncation_does_not_break_statistics() {
    // Use a sum_cutoff large enough that low-mass branches are dropped
    // and the remainder is renormalized. Statistics should still be
    // close to the pure-tableau ground truth.
    let shots = 8000;
    let sum_cutoff = 1e-4;
    let sum = run_sum(3, shots, sum_cutoff, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        for q in 0..3 {
            t.depolarize(q, 0.05);
            t.loss_channel(q, 0.05);
        }
    });
    let pure = run_pure(3, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        for q in 0..3 {
            t.depolarize(q, 0.05);
            t.loss_channel(q, 0.05);
        }
    });
    // Slightly looser tolerance to allow for truncation-induced bias.
    assert_distributions_match(&sum, &pure, 0.1, "truncation_does_not_break_statistics");
}

// ---------------------------------------------------------------------------
// Mid-circuit measurement
// ---------------------------------------------------------------------------
//
// `GeneralizedTableauSum::measure(addr0)` branches the sum into the possible
// outcomes (and a separate lost-qubit bucket) rather than collapsing to a
// single trajectory. After a mid-circuit measurement, sampling the resulting
// sum + `measure_all` on the chosen branch must reproduce the same joint
// distribution as the pure backend, which collapses stochastically per shot
// via `LossyMeasure::measure`.

#[test]
fn mid_circuit_measure_on_zero_state_is_deterministic() {
    // |0⟩ is a Z eigenstate (case b, no branching). Both backends must
    // return Some(false) for every shot.
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        let _ = t.measure(0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(false)));
    assert!(pure.iter().all(|s| s[0] == Some(false)));
}

#[test]
fn mid_circuit_measure_after_x_is_deterministic_one() {
    // X|0⟩ = |1⟩ — also a Z eigenstate (case b), deterministic outcome 1.
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.x(0);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.x(0);
        let _ = t.measure(0);
    });
    assert!(sum.iter().all(|s| s[0] == Some(true)));
    assert!(pure.iter().all(|s| s[0] == Some(true)));
}

#[test]
fn mid_circuit_measure_after_h_is_unbiased() {
    // H|0⟩ = |+⟩ is not a Z eigenstate → case a, sum branches 50/50.
    // Sampled marginal must match pure backend's 50/50.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        let _ = t.measure(0);
    });
    assert_distributions_match(&sum, &pure, 0.04, "mid_circuit_measure_after_h_is_unbiased");
}

#[test]
fn mid_circuit_measure_then_h_back_is_unbiased() {
    // After the measurement the qubit is in |0⟩ or |1⟩ (each 50%);
    // H sends both to |±⟩ which measure to 50/50 in Z. Total: 50/50.
    // Verifies that post-measurement state evolution is consistent
    // across both backends (gates must apply correctly to the branched sum).
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.measure(0);
        t.h(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        let _ = t.measure(0);
        t.h(0);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_then_h_back_is_unbiased",
    );
}

#[test]
fn mid_circuit_measure_bell_pair_correlates_outcomes() {
    // Bell pair (|00⟩+|11⟩)/√2, then measure q0 mid-circuit. The post-state
    // collapses both qubits to the same basis vector, so final measure_all
    // yields only (0,0) and (1,1) (each ~50%).
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.measure(0);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        let _ = t.measure(0);
    });
    assert!(
        sum.iter().all(|s| s[0] == s[1]),
        "sum: Bell pair must remain correlated through mid-circuit measure"
    );
    assert!(
        pure.iter().all(|s| s[0] == s[1]),
        "pure: Bell pair must remain correlated through mid-circuit measure"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_bell_pair_correlates_outcomes",
    );
}

#[test]
fn mid_circuit_measure_bell_pair_then_h_q1_decorrelates() {
    // After measuring q0, q1 is in a definite state (|0⟩ or |1⟩, 50/50).
    // H on q1 puts it in |+⟩ or |−⟩ — measure_all gives independent 50/50.
    // Joint: q0 ∈ {0,1} 50/50, q1 ∈ {0,1} 50/50, independent → 4 bins ×25%.
    let shots = 8000;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.measure(0);
        t.h(1);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        let _ = t.measure(0);
        t.h(1);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.05,
        "mid_circuit_measure_bell_pair_then_h_q1_decorrelates",
    );
}

#[test]
fn mid_circuit_measure_consecutive_is_idempotent() {
    // Two measurements on the same qubit with nothing between them must
    // agree (case a then case b). The second is a no-op on the post-state.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.measure(0);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        let _ = t.measure(0);
        let _ = t.measure(0);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_consecutive_is_idempotent",
    );
}

#[test]
fn mid_circuit_measure_on_lost_qubit_stays_lost() {
    // Loss at p=1.0 marks q0 lost; subsequent measure must leave the
    // qubit lost (`None`) on both backends.
    let shots = 1000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.loss_channel(0, 1.0);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.loss_channel(0, 1.0);
        let _ = t.measure(0);
    });
    assert!(sum.iter().all(|s| s[0].is_none()));
    assert!(pure.iter().all(|s| s[0].is_none()));
}

// ---------------------------------------------------------------------------
// Mid-circuit measurement with multiple incoming branches
//
// Noise channels create multiple sum entries; measure() must branch each one
// independently and keep the joint distribution correct after merge.
// ---------------------------------------------------------------------------

#[test]
fn mid_circuit_measure_after_depolarize_on_plus_state() {
    // h(0) puts q0 in |+⟩; depolarize creates 4 branches (I/X/Y/Z applied),
    // then measure on the resulting mixture. |+⟩ and X|+⟩=|+⟩ measure 50/50,
    // Y|+⟩=-i|−⟩ and Z|+⟩=|−⟩ also measure 50/50, so total is 50/50 — but
    // the test value is in covering the multi-branch-into-measure path.
    let shots = 8000;
    let p = 0.3_f64;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.depolarize(0, p);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.depolarize(0, p);
        let _ = t.measure(0);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_after_depolarize_on_plus_state",
    );
}

#[test]
fn mid_circuit_measure_bell_pair_after_depolarize() {
    // Bell pair, depolarize on q0 (4 branches), then measure(q0) mid-circuit,
    // then measure_all. The depolarize breaks the Bell correlation: with
    // probability 2p/3 the two qubits disagree, otherwise they agree. The
    // mid-circuit measure on q0 must propagate that disagreement correctly
    // through the rest of the circuit.
    let shots = 8000;
    let p = 0.25_f64;
    let sum = run_sum(2, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.depolarize(0, p);
        t.measure(0);
    });
    let pure = run_pure(2, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.depolarize(0, p);
        let _ = t.measure(0);
    });
    // Sanity: P(q0 ≠ q1) should be ~2p/3 on both backends.
    let disagree_sum = sum.iter().filter(|s| s[0] != s[1]).count() as f64 / shots as f64;
    let disagree_pure = pure.iter().filter(|s| s[0] != s[1]).count() as f64 / shots as f64;
    let expected = 2.0 * p / 3.0;
    assert!(
        (disagree_sum - expected).abs() < 0.04,
        "sum P(disagree)={disagree_sum:.4}, expected {expected:.4}"
    );
    assert!(
        (disagree_pure - expected).abs() < 0.04,
        "pure P(disagree)={disagree_pure:.4}, expected {expected:.4}"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.05,
        "mid_circuit_measure_bell_pair_after_depolarize",
    );
}

#[test]
fn mid_circuit_measure_after_loss_includes_none_bucket() {
    // Loss at p=0.4 splits the sum into present/lost branches; measure
    // operates on both. The lost branch passes through as None; the
    // present branch (originally |+⟩) gives 50/50.
    // Expected per qubit: P(None) = 0.4, P(0) = 0.3, P(1) = 0.3.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.loss_channel(0, 0.4);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        t.loss_channel(0, 0.4);
        let _ = t.measure(0);
    });
    let lost_sum = sum.iter().filter(|s| s[0].is_none()).count() as f64 / shots as f64;
    let lost_pure = pure.iter().filter(|s| s[0].is_none()).count() as f64 / shots as f64;
    assert!(
        (lost_sum - 0.4).abs() < 0.04,
        "sum P(lost)={lost_sum:.4}, expected 0.4"
    );
    assert!(
        (lost_pure - 0.4).abs() < 0.04,
        "pure P(lost)={lost_pure:.4}, expected 0.4"
    );
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_after_loss_includes_none_bucket",
    );
}

#[test]
fn mid_circuit_measure_then_depolarize_then_measure() {
    // Interleave measure with noise: first measure collapses to |0⟩/|1⟩,
    // depolarize creates branches, second measure resolves each. Stresses
    // the storage layer's ability to handle measure → re-bucket → noise
    // → re-bucket cycles.
    let shots = 8000;
    let sum = run_sum(1, shots, 1e-12, |t| {
        t.h(0);
        t.measure(0);
        t.depolarize(0, 0.2);
        t.measure(0);
    });
    let pure = run_pure(1, shots, |t| {
        t.h(0);
        let _ = t.measure(0);
        t.depolarize(0, 0.2);
        let _ = t.measure(0);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.04,
        "mid_circuit_measure_then_depolarize_then_measure",
    );
}

#[test]
fn mid_circuit_measure_three_qubit_with_noise() {
    // Larger circuit: GHZ-like state, mid-circuit measure on the middle
    // qubit, then more gates and a final measure_all. Verifies that the
    // sum's bookkeeping (fingerprint cache, merge, normalize) survives a
    // realistic measurement-aware noisy circuit.
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        t.depolarize(1, 0.1);
        t.measure(1);
        t.h(2);
        t.depolarize(0, 0.1);
    });
    let pure = run_pure(3, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(1, 2);
        t.depolarize(1, 0.1);
        let _ = t.measure(1);
        t.h(2);
        t.depolarize(0, 0.1);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.08,
        "mid_circuit_measure_three_qubit_with_noise",
    );
}

#[test]
fn mid_circuit_measure_case_a_many_entries_with_other_outcome_retained() {
    // Regression test for the case-a path of GeneralizedTableauSum::measure
    // when (a) the sum holds many entries and (b) both projection outcomes
    // are above sum_cutoff on every iteration, so both project_case_a calls
    // fire for every entry.
    //
    // This stresses the shared-across-iterations scratch refactor:
    //   - `scratch.coeff_map` must be cleared between iterations, else stale
    //     entries from the previous tableau leak into compute_overlap_case_a
    //     and the partition inside project_case_a.
    //   - `scratch_other_outcome.coeff_map` must be repopulated (via
    //     clone_from) before the other-outcome project_case_a — otherwise
    //     the branch keeps its forked, un-projected coefficients.
    //
    // GHZ + repeated depolarize on q1/q2 keeps q0 in superposition for every
    // sum entry (so measure(0) is case-a everywhere) while multiplying the
    // entry count. Post-measurement gates probe the projected state on both
    // branches; a wrong projection on the other-outcome branch shows up as
    // statistical divergence from the per-shot pure backend.
    let shots = 8000;
    let sum = run_sum(3, shots, 1e-12, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(0, 2);
        for _ in 0..3 {
            t.depolarize(1, 0.18);
            t.depolarize(2, 0.18);
        }
        t.measure(0);
        t.h(1);
        t.h(2);
    });
    let pure = run_pure(3, shots, |t| {
        t.h(0);
        t.cnot(0, 1);
        t.cnot(0, 2);
        for _ in 0..3 {
            t.depolarize(1, 0.18);
            t.depolarize(2, 0.18);
        }
        let _ = t.measure(0);
        t.h(1);
        t.h(2);
    });
    assert_distributions_match(
        &sum,
        &pure,
        0.08,
        "mid_circuit_measure_case_a_many_entries_with_other_outcome_retained",
    );
}
