use std::sync::Mutex;

use ppvm_runtime::prelude::ACMap;
use ppvm_runtime::prelude::PauliStorage;
use ppvm_runtime::traits::{Coefficient, PauliWordTrait, Strategy};

/// Truncation strategy that combines a coefficient threshold with a hard cap on `|P|`.
///
/// First removes entries below `min_threshold` (identical to `CoefficientThreshold`).
/// Then, if the map still has more than `target` entries, prunes the excess.
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
    /// Entries with `|coeff| < min_threshold` are always removed.
    pub min_threshold: f64,
}

impl Default for Budget {
    fn default() -> Self {
        Budget { target: usize::MAX, min_threshold: 1e-12 }
    }
}

impl Strategy for Budget {
    fn capacity(&self, n_qubits: usize) -> usize {
        // Use target as capacity hint when it's a real bound; fall back to the
        // CoefficientThreshold heuristic when target is effectively unlimited.
        if self.target < usize::MAX / 2 {
            self.target
        } else {
            n_qubits * 10
        }
    }

    fn truncate<S, V, H, M, W>(&self, map: &mut M)
    where
        S: PauliStorage,
        V: Coefficient,
        H: std::hash::BuildHasher + Clone + Default,
        W: PauliWordTrait,
        M: ACMap<S, V, H, W>,
    {
        // Step 1: threshold — identical to CoefficientThreshold(min_threshold).
        map.retain(|_, v| !v.cutoff(self.min_threshold));

        if map.len() <= self.target {
            return;
        }

        // Step 2: find the magnitude threshold t* such that at most `target`
        // entries have |v| >= t*, then retain those entries.
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
        let mut hi = self.min_threshold.max(f64::MIN_POSITIVE) * 2.0;
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
    use ppvm_runtime::config::fxhash::ByteF64;
    use ppvm_runtime::prelude::{PauliSum, PauliWord, PhasedPauliWord};
    use ppvm_runtime::strategy::CoefficientThreshold;

    use crate::lindblad::{CollapseOp, JumpOp, LindbladOp, RateMatrix};
    use crate::solve::{SolverConfig, solve};

    type W1 = PauliWord<[u8; 1], fxhash::FxBuildHasher>;
    type BudgetConfig = ByteF64<1, Budget>;
    type ThreshConfig = ByteF64<1, CoefficientThreshold>;

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

    /// Budget must keep |P| ≤ target at every save point.
    #[test]
    fn budget_limits_size() {
        // n=1, target=1: after every truncation step, at most 1 entry survives.
        let strat = Budget { target: 1, min_threshold: 0.0 };
        let initial: PauliSum<BudgetConfig> =
            PauliSum::builder().n_qubits(1).strategy(strat).build();
        // Start with a non-trivial state by using solve with initial P = Z.
        let mut initial_z: PauliSum<BudgetConfig> =
            PauliSum::builder().n_qubits(1).strategy(strat).build();
        initial_z += ("Z", 1.0_f64);
        let _ = initial; // ensure default builds

        let lop = lindblad_lowering();
        let save_at: Vec<f64> = (1..=5).map(|i| i as f64 * 0.05).collect();
        let config = SolverConfig::default();

        // solve returns saved states; we measure |P| at each save.
        let (_, sizes) = solve(
            None, &lop, &initial_z, (0.0, 0.25), &save_at,
            |_, p| p.data().len(),
            config,
        );

        for (i, &sz) in sizes.iter().enumerate() {
            assert!(sz <= 1, "save point {i}: |P| = {sz} exceeds target=1");
        }
    }

    /// When |P| ≤ target, Budget must give the same trajectory as CoefficientThreshold.
    #[test]
    fn budget_matches_threshold() {
        // n=1, spontaneous emission. |P| stays ≤ 4 << target=100.
        // Both strategies use the same min_threshold, so results must be identical.
        let threshold = 1e-8_f64;
        let save_at: Vec<f64> = vec![0.1, 0.2, 0.3];
        let config = SolverConfig { rtol: threshold, ..SolverConfig::default() };

        // CoefficientThreshold reference.
        let mut initial_ct: PauliSum<ThreshConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(CoefficientThreshold(threshold))
            .build();
        initial_ct += ("Z", 1.0_f64);
        let (_, ct_z) = solve(
            None, &lindblad_lowering_thresh(), &initial_ct, (0.0, 0.3), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            config,
        );

        // Budget with the same threshold but a generous target.
        let mut initial_b: PauliSum<BudgetConfig> = PauliSum::builder()
            .n_qubits(1)
            .strategy(Budget { target: 100, min_threshold: threshold })
            .build();
        initial_b += ("Z", 1.0_f64);
        let config2 = SolverConfig { rtol: threshold, ..SolverConfig::default() };
        let (_, b_z) = solve(
            None, &lindblad_lowering(), &initial_b, (0.0, 0.3), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            config2,
        );

        assert_eq!(ct_z.len(), b_z.len());
        for (i, (&ct, &b)) in ct_z.iter().zip(b_z.iter()).enumerate() {
            assert!(
                (ct - b).abs() < 1e-14,
                "save {i}: CoefficientThreshold={ct}, Budget={b}, diff={:.2e}",
                (ct - b).abs()
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
            .strategy(Budget { target: 200, min_threshold: 1e-12 })
            .build();
        initial_b += ("Z", 1.0_f64);
        let config2 = SolverConfig::default();
        let (_, bud_z): (_, Vec<f64>) = solve(
            None, &lindblad_lowering(), &initial_b, (0.0, 0.5), &save_at,
            |_, p| { use ppvm_runtime::prelude::Trace; p.data().trace(&W1::from("Z")) },
            config2,
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
