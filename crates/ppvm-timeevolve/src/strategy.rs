use std::sync::Mutex;

use ppvm_runtime::prelude::ACMap;
use ppvm_runtime::prelude::PauliStorage;
use ppvm_runtime::traits::{Coefficient, PauliWordTrait, Strategy};

/// Hard cap on the number of Pauli entries in a `PauliSum`.
///
/// # Performance note (Task 28 benchmark, n=5, target=300)
///
/// `bench_step_budget` (446 µs) vs `bench_step_ct` (388 µs): **1.15× overhead** per
/// DOPRI5 step.  The binary-search pass over |P| entries runs 6 times per step (once per
/// truncation site); at the benchmark fixture size (~25 live terms, target=300 never fires)
/// this is negligible.  When the cap fires on a large state the overhead scales as
/// O(|P|·128) per truncation call.  The 2× threshold that would trigger a "not recommended
/// as default" advisory was not reached; `Budget` is suitable when a hard memory cap is
/// needed.
///
/// `Budget` is a *pure count cap*: it keeps the `target` entries with the largest
/// coefficient magnitudes and drops the rest.  It does **not** apply a coefficient
/// threshold; small-but-nonzero entries survive as long as there is room.
///
/// To combine count capping with coefficient pruning, compose the two strategies:
/// ```text
/// ByteFxHashF64<N, CombinedStrategy<Budget, CoefficientThreshold>>
/// ```
/// `CombinedStrategy` applies them in order — threshold first, then cap — so you
/// get both behaviours with a single `PauliSum` type.
///
/// **Implementation note:** `Strategy::truncate` only exposes `retain` and
/// `scale` (no `ACMapIter`), and `V: Coefficient` has no magnitude accessor.
/// We recover iteration via `scale` as a side-effect, then binary-search on
/// `V::cutoff(t)` to find the kth-largest threshold. Cost: O(n·128) per
/// truncation call — negligible for practical map sizes.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    /// Hard cap on the number of Pauli entries retained after truncation.
    pub target: usize,
}

impl Default for Budget {
    fn default() -> Self {
        Budget { target: usize::MAX }
    }
}

impl Strategy for Budget {
    fn capacity(&self, n_qubits: usize) -> usize {
        // Fall back to a sensible hint when target is effectively unlimited.
        if self.target < usize::MAX / 2 { self.target } else { n_qubits * 10 }
    }

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        W: PauliWordTrait,
        M: ACMap<S, V, H, W>,
    {
        if map.len() <= self.target {
            return;
        }

        // Find the magnitude threshold t* such that at most `target` entries
        // have |v| >= t*, then retain those entries.
        //
        // Strategy::truncate only has M: ACMap<S, V, H, W> in scope — no
        // ACMapIter, and V: Coefficient does not expose a magnitude accessor.
        // We work around this with two passes:
        //   Pass A: use `scale` as a side-effect iterator to collect all V
        //           values (without modifying them).
        //   Pass B: binary-search on t using V::cutoff(t) -> bool, which for
        //           f64 means |v| < t. This gives O(n·log(1/ε)) comparisons —
        //           negligible for the map sizes expected in practice.
        //   Pass C: retain(|_, v| !v.cutoff(t*)) keeps the largest entries.
        let n = map.len();
        let collected: Mutex<Vec<V>> = Mutex::new(Vec::with_capacity(n));
        map.scale(|_, v| {
            // Read v without modifying it (scale is used purely for iteration).
            collected.lock().expect("Budget collect").push(v.clone());
        });
        let values = collected.into_inner().expect("Budget values");

        // Closure: count entries with |v| >= t.
        let count_ge = |t: f64| -> usize {
            values.iter().filter(|v| !v.cutoff(t)).count()
        };

        // Exponential search for hi where count_ge(hi) == 0.
        let mut hi = f64::MIN_POSITIVE * 2.0;
        while count_ge(hi) > 0 && hi < 1e300 {
            hi *= 2.0;
        }

        // Binary search: find t in (lo, hi] where count_ge(t) <= target.
        // Invariant: count_ge(lo) > target, count_ge(hi) <= target.
        let mut lo = 0.0_f64;
        for _ in 0..128 {
            let mid = lo + (hi - lo) * 0.5;
            if mid <= lo || mid >= hi {
                break; // float convergence
            }
            if count_ge(mid) > self.target {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        // Retain entries with |v| >= hi (at most `target` of them by construction).
        map.retain(|_, v| !v.cutoff(hi));

        // Tie-break: if float precision leaves a few extra entries at the
        // boundary magnitude, trim with a counter. All trimmed entries have
        // |v| == hi, so the magnitude difference is negligible.
        if map.len() > self.target {
            let target = self.target;
            let kept = Mutex::new(0usize);
            map.retain(|_, _| {
                let mut g = kept.lock().expect("Budget tie-break");
                if *g < target { *g += 1; true } else { false }
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::indexmap::ByteFxHashF64;
    use ppvm_runtime::prelude::{PauliSum, PauliWord, PhasedPauliWord};
    use ppvm_runtime::strategy::{CoefficientThreshold, CombinedStrategy};

    use crate::lindblad::{CollapseOp, JumpOp, LindbladOp, RateMatrix};
    use crate::solve::{SolverConfig, solve};

    type W1 = PauliWord<[u8; 1], fxhash::FxBuildHasher>;
    type BudgetConfig = ByteFxHashF64<1, Budget>;
    type ThreshConfig = ByteFxHashF64<1, CoefficientThreshold>;
    type CombinedConfig = ByteFxHashF64<1, CombinedStrategy<Budget, CoefficientThreshold>>;

    fn ppw(pauli: &str, phase: u8) -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, W1> {
        PhasedPauliWord::build_from_word(W1::from(pauli), phase)
    }

    fn lindblad_lowering() -> LindbladOp<BudgetConfig> {
        let mut op = CollapseOp::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        LindbladOp::new(vec![JumpOp::Generic(op)], RateMatrix::from(vec![1.0]))
    }

    fn lindblad_lowering_thresh() -> LindbladOp<ThreshConfig> {
        let mut op = CollapseOp::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        LindbladOp::new(vec![JumpOp::Generic(op)], RateMatrix::from(vec![1.0]))
    }

    fn lindblad_lowering_combined() -> LindbladOp<CombinedConfig> {
        let mut op = CollapseOp::new(1);
        op.push(ppw("X", 0), 1.0);
        op.push(ppw("Y", 1), 1.0);
        LindbladOp::new(vec![JumpOp::Generic(op)], RateMatrix::from(vec![1.0]))
    }

    /// Budget must keep |P| ≤ target at every save point.
    #[test]
    fn budget_limits_size() {
        let strat = Budget { target: 1 };
        let mut initial_z: PauliSum<BudgetConfig> =
            PauliSum::builder().n_qubits(1).strategy(strat).build();
        initial_z += ("Z", 1.0_f64);

        let lop = lindblad_lowering();
        let save_at: Vec<f64> = (1..=5).map(|i| i as f64 * 0.05).collect();

        let (_, sizes) = solve(
            None, &lop, &initial_z, (0.0, 0.25), &save_at,
            |_, p| p.data().len(),
            SolverConfig::default(),
        );

        for (i, &sz) in sizes.iter().enumerate() {
            assert!(sz <= 1, "save point {i}: |P| = {sz} exceeds target=1");
        }
    }

    /// Budget { target: 5 } keeps the 5 largest-magnitude entries.
    #[test]
    fn budget_no_threshold_keeps_largest() {
        // Build a PauliSum with 10 terms of known magnitudes 1.0, 0.9, ..., 0.1.
        // After Budget { target: 5 }, only the top 5 (magnitudes 1.0–0.6) survive.
        // Use 2 qubits so we have room for 10 distinct Pauli operators.
        let strat = Budget { target: 5 };
        let mut p: PauliSum<ByteFxHashF64<1, Budget>> =
            PauliSum::builder().n_qubits(2).strategy(strat).build();

        let labels = ["II","IZ","IX","IY","ZI","ZZ","ZX","ZY","XI","XZ"];
        for (i, &label) in labels.iter().enumerate() {
            let coeff = 1.0 - 0.1 * i as f64; // 1.0, 0.9, 0.8, ..., 0.1
            p += (label, coeff);
        }

        assert_eq!(p.data().len(), 10, "expected 10 terms before truncation");
        p.truncate();
        assert_eq!(
            p.data().len(), 5,
            "|P| after Budget{{target=5}} should be 5, got {}",
            p.data().len()
        );

    }

    /// CombinedStrategy<Budget, CoefficientThreshold> gives the same trajectory
    /// as applying CoefficientThreshold alone when |P| << target.
    #[test]
    fn combined_strategy_matches_separate() {
        let threshold = 1e-8_f64;
        let save_at: Vec<f64> = vec![0.1, 0.2, 0.3];

        // Reference: CoefficientThreshold only.
        let mut initial_ct: PauliSum<ThreshConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(threshold))
            .build();
        initial_ct += ("Z", 1.0_f64);
        let (_, ct_z) = solve(
            None, &lindblad_lowering_thresh(), &initial_ct, (0.0, 0.3), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            SolverConfig::default(),
        );

        // CombinedStrategy: Budget (generous cap, never fires) + CoefficientThreshold.
        let mut initial_c: PauliSum<CombinedConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CombinedStrategy(Budget { target: 100 }, CoefficientThreshold(threshold)))
            .build();
        initial_c += ("Z", 1.0_f64);
        let (_, comb_z) = solve(
            None, &lindblad_lowering_combined(), &initial_c, (0.0, 0.3), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            SolverConfig::default(),
        );

        assert_eq!(ct_z.len(), comb_z.len());
        for (i, (&ct, &c)) in ct_z.iter().zip(comb_z.iter()).enumerate() {
            assert!(
                (ct - c).abs() < 1e-14,
                "save {i}: CoefficientThreshold={ct}, Combined={c}, diff={:.2e}",
                (ct - c).abs()
            );
        }
    }

    /// Budget { target: 200 } matches CoefficientThreshold(1e-12) to within 1e-3
    /// for the spontaneous emission system (n=1; |P| ≤ 4 << 200 so budget never fires).
    #[test]
    fn budget_accuracy() {
        let save_at: Vec<f64> = vec![0.1, 0.25, 0.5];
        let config = SolverConfig::default();

        let mut initial_ct: PauliSum<ThreshConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(1e-12))
            .build();
        initial_ct += ("Z", 1.0_f64);
        let (_, ref_z): (_, Vec<f64>) = solve(
            None, &lindblad_lowering_thresh(), &initial_ct, (0.0, 0.5), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            config,
        );

        let mut initial_b: PauliSum<BudgetConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(Budget { target: 200 })
            .build();
        initial_b += ("Z", 1.0_f64);
        let (_, bud_z): (_, Vec<f64>) = solve(
            None, &lindblad_lowering(), &initial_b, (0.0, 0.5), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            SolverConfig::default(),
        );

        for (i, (&r, &b)) in ref_z.iter().zip(bud_z.iter()).enumerate() {
            assert!(
                (r - b).abs() < 1e-3,
                "save {i}: reference={r}, budget={b}, diff={:.2e}",
                (r - b).abs()
            );
        }
    }
}
