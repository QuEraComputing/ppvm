//! Statistical equivalence tests between the noise-aware `Sampler` (sum
//! backend) and the pure `GeneralizedTableau` simulator.
//!
//! The sum backend evolves a probability-weighted collection of branches
//! through `loss_channel` / `depolarize` (one branch per error outcome)
//! and then samples a branch + measurement pair per shot. The pure
//! tableau applies each channel stochastically inside a single trajectory.
//! In the limit of many shots both must yield the same joint distribution
//! over `Vec<Option<bool>>` measurement outcomes.
//!
//! Each test runs N shots on both backends with deterministic seeds and
//! checks that the total variation distance between the empirical
//! distributions is below a finite-sample threshold.

use std::collections::{HashMap, HashSet};

use ppvm_runtime::config::fxhash::ByteF64;
use ppvm_runtime::traits::{Clifford, CliffordExtensions, Depolarizing, LossChannel};
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
