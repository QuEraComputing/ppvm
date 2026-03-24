use std::ops::{Mul, MulAssign};

use ppvm_runtime::prelude::{
    ACMapAddAssign, ACMapBase, ACMapIter, Config, PauliSum, PauliWord, PauliWordTrait,
    PhasedPauliWord, Trace,
};

use crate::lindblad::{LindbladOp, rhs, rhs_into};
use crate::solve::{SolverCache, SolverConfig};

// ---------------------------------------------------------------------------
// Dormand-Prince 4(5) Butcher tableau
// Reference: Hairer et al. "Solving ODEs I", Table 5.2
// ---------------------------------------------------------------------------

// Stage coefficients a_{i,j} (1-indexed, lower-triangular)
const A21: f64 = 1.0 / 5.0;
const A31: f64 = 3.0 / 40.0;
const A32: f64 = 9.0 / 40.0;
const A41: f64 = 44.0 / 45.0;
const A42: f64 = -56.0 / 15.0;
const A43: f64 = 32.0 / 9.0;
const A51: f64 = 19372.0 / 6561.0;
const A52: f64 = -25360.0 / 2187.0;
const A53: f64 = 64448.0 / 6561.0;
const A54: f64 = -212.0 / 729.0;
const A61: f64 = 9017.0 / 3168.0;
const A62: f64 = -355.0 / 33.0;
const A63: f64 = 46732.0 / 5247.0;
const A64: f64 = 49.0 / 176.0;
const A65: f64 = -5103.0 / 18656.0;

// 5th-order b weights (= a_{7,*} via FSAL; b_7 = 0 so k7 does not appear in y_new)
const B1: f64 = 35.0 / 384.0;
// B2 = 0
const B3: f64 = 500.0 / 1113.0;
const B4: f64 = 125.0 / 192.0;
const B5: f64 = -2187.0 / 6784.0;
const B6: f64 = 11.0 / 84.0;

// Error coefficients e_i = b_i - b*_i (7 entries, e_2 = 0)
// b*: 4th-order embedded solution weights
// b*_1=5179/57600, b*_3=7571/16695, b*_4=393/640,
// b*_5=-92097/339200, b*_6=187/2100, b*_7=1/40
const E1: f64 = 71.0 / 57600.0;
// E2 = 0
const E3: f64 = -71.0 / 16695.0;
const E4: f64 = 71.0 / 1920.0;
const E5: f64 = -17253.0 / 339200.0;
const E6: f64 = 22.0 / 525.0;
const E7: f64 = -1.0 / 40.0;

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

pub(crate) enum StepResult {
    Accept { h_new: f64 },
    Reject { h_new: f64 },
}

// ---------------------------------------------------------------------------
// Helper: accumulate `target += scale * source` term by term
// ---------------------------------------------------------------------------

fn add_scaled<T: Config>(target: &mut PauliSum<T>, source: &PauliSum<T>, scale: f64)
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: Copy + std::ops::AddAssign + std::ops::Mul<Output = T::Coeff>,
    T::PauliWordType: Clone,
    f64: Into<T::Coeff>,
{
    let s: T::Coeff = scale.into();
    for (w, c) in source.data().iter() {
        // PhasedPauliWord::mul_assign updates xbits/zbits but not hash_cache.
        // Rehash here so that words from k-vector multiplication entries hash
        // consistently with words constructed from strings, preventing duplicate
        // map entries for the same logical Pauli key.
        let mut key = w.clone();
        key.rehash();
        *target += (key, s * *c);
    }
}

// ---------------------------------------------------------------------------
// h0 auto-estimation (Hairer et al. "Solving ODEs I", §II.4)
// ---------------------------------------------------------------------------

/// Estimate the initial step size for the ODE solver.
///
/// If `config.h0` is set, returns that value clamped to `[hmin, hmax_eff]`.
/// Otherwise applies the 5-step Hairer procedure and returns a suitable h0.
/// `hmax_eff = config.hmax.min(t_span.1 - t_span.0)`.
#[allow(dead_code)] // called from solve.rs (Task 9)
pub(crate) fn estimate_h0<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    y0: &PauliSum<T>,
    t_span: (f64, f64),
    config: &SolverConfig,
) -> f64
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + ACMapBase
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign
        + Copy
        + std::ops::Mul<Output = T::Coeff>
        + std::iter::Sum
        + Into<f64>
        + Send,
    T::PauliWordType: Clone
        + std::borrow::Borrow<PauliWord<T::Storage, T::BuildHasher>>
        + Send
        + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
    for<'a> T::Map: Trace<'a, T::PauliWordType, Output = T::Coeff>,
{
    let hmax_eff = config.hmax.min(t_span.1 - t_span.0);

    if let Some(h0) = config.h0 {
        return h0.clamp(config.hmin, hmax_eff);
    }

    // Step 1: norms of y0 and f0
    let d0 = y0.overlap(y0).into().sqrt();
    let f0 = rhs(ham, lindblad, y0);
    let d1 = f0.overlap(&f0).into().sqrt();

    // Step 2: rough initial h0 (fallback to 1e-6 if scales are too small)
    let h0 = if d0 < 1e-5 || d1 < 1e-5 { 1e-6 } else { 0.01 * d0 / d1 };

    // Step 3: one explicit Euler step, compute derivative there
    let mut y1 = y0.clone();
    add_scaled(&mut y1, &f0, h0);
    let f1 = rhs(ham, lindblad, &y1);

    // Step 4: h1 from second-derivative estimate
    // d2 = ||f1 - f0|| / h0
    let mut df = f1;
    add_scaled(&mut df, &f0, -1.0);
    let d2 = df.overlap(&df).into().sqrt() / h0;

    let h1 = if d1 <= 1e-5 && d2 <= 1e-5 {
        (1e-6_f64).max(h0 * 1e-3)
    } else {
        (0.01 / d1.max(d2)).powf(0.2)
    };

    // Step 5: pick the smallest of the three candidates, clamped to [hmin, hmax_eff]
    let h_init = (100.0 * h0).min(h1).min(hmax_eff);
    h_init.clamp(config.hmin, hmax_eff)
}

// ---------------------------------------------------------------------------
// Single adaptive Dormand-Prince step
// ---------------------------------------------------------------------------

/// Attempt one DOPRI5 step from `y` with step size `dt`.
///
/// On entry `cache.k[0]` must hold `rhs(y)` (FSAL carry-over from the previous step).
/// On `Accept`, `cache.y_scratch` holds the 5th-order solution and `cache.k[6]` holds
/// `rhs(y_scratch)` (the next k1).  The caller is responsible for swapping state and
/// advancing the FSAL index (`cache.k.swap(0, 6)`).
/// Returns `Accept` if the local error is within tolerance, `Reject` otherwise.
pub(crate) fn step<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    y: &PauliSum<T>,
    dt: f64,
    config: &SolverConfig,
    cache: &mut SolverCache<T>,
) -> StepResult
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>
        + ACMapBase
        + Clone
        + Send
        + Sync,
    T::Coeff: std::ops::AddAssign
        + Copy
        + std::ops::Mul<Output = T::Coeff>
        + std::iter::Sum
        + Into<f64>
        + Send,
    T::PauliWordType: Clone
        + std::borrow::Borrow<PauliWord<T::Storage, T::BuildHasher>>
        + Send
        + Sync,
    T::BuildHasher: Send + Sync,
    T::Strategy: Send + Sync,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
    for<'a> T::Map: Trace<'a, T::PauliWordType, Output = T::Coeff>,
{
    // Destructure cache fields so the borrow checker can track independent borrows.
    let SolverCache { k, y_scratch, err } = cache;

    // Stage 2: y_scratch = y + dt*A21*k[0];  k[1] = rhs(y_scratch)
    {
        let (lo, hi) = k.split_at_mut(1);
        y_scratch.data_mut().clone_from(y.data());
        add_scaled(y_scratch, &lo[0], dt * A21);
        rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
    }

    // Stage 3: y_scratch = y + dt*(A31*k[0] + A32*k[1]);  k[2] = rhs(y_scratch)
    {
        let (lo, hi) = k.split_at_mut(2);
        y_scratch.data_mut().clone_from(y.data());
        add_scaled(y_scratch, &lo[0], dt * A31);
        add_scaled(y_scratch, &lo[1], dt * A32);
        rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
    }

    // Stage 4
    {
        let (lo, hi) = k.split_at_mut(3);
        y_scratch.data_mut().clone_from(y.data());
        add_scaled(y_scratch, &lo[0], dt * A41);
        add_scaled(y_scratch, &lo[1], dt * A42);
        add_scaled(y_scratch, &lo[2], dt * A43);
        rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
    }

    // Stage 5
    {
        let (lo, hi) = k.split_at_mut(4);
        y_scratch.data_mut().clone_from(y.data());
        add_scaled(y_scratch, &lo[0], dt * A51);
        add_scaled(y_scratch, &lo[1], dt * A52);
        add_scaled(y_scratch, &lo[2], dt * A53);
        add_scaled(y_scratch, &lo[3], dt * A54);
        rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
    }

    // Stage 6
    {
        let (lo, hi) = k.split_at_mut(5);
        y_scratch.data_mut().clone_from(y.data());
        add_scaled(y_scratch, &lo[0], dt * A61);
        add_scaled(y_scratch, &lo[1], dt * A62);
        add_scaled(y_scratch, &lo[2], dt * A63);
        add_scaled(y_scratch, &lo[3], dt * A64);
        add_scaled(y_scratch, &lo[4], dt * A65);
        rhs_into(ham, lindblad, y_scratch, &mut hi[0]);
    }

    // 5th-order update: y_scratch = y + dt*(b1*k[0] + b3*k[2] + b4*k[3] + b5*k[4] + b6*k[5])
    // b_2 = 0 and b_7 = 0 (FSAL: k7 not used in y_new)
    y_scratch.data_mut().clone_from(y.data());
    add_scaled(y_scratch, &k[0], dt * B1);
    add_scaled(y_scratch, &k[2], dt * B3);
    add_scaled(y_scratch, &k[3], dt * B4);
    add_scaled(y_scratch, &k[4], dt * B5);
    add_scaled(y_scratch, &k[5], dt * B6);

    // Truncate the new state before FSAL so sub-threshold Pauli strings do not
    // accumulate across steps. Doing this here (before Stage 7) keeps k[6]
    // consistent with the truncated state used at the start of the next step.
    y_scratch.truncate();

    // Stage 7 (FSAL): k[6] = rhs(y_scratch)
    rhs_into(ham, lindblad, y_scratch, &mut k[6]);

    // Error estimate: e = dt * (e1*k[0] + e3*k[2] + e4*k[3] + e5*k[4] + e6*k[5] + e7*k[6])
    // e_2 = 0
    err.data_mut().clear();
    add_scaled(err, &k[0], dt * E1);
    add_scaled(err, &k[2], dt * E3);
    add_scaled(err, &k[3], dt * E4);
    add_scaled(err, &k[4], dt * E5);
    add_scaled(err, &k[5], dt * E6);
    add_scaled(err, &k[6], dt * E7);

    // Error norm: err = ||e|| / (atol + rtol * ||y||)
    let err_norm_sq: f64 = err.overlap(err).into();
    let y_norm_sq: f64 = y.overlap(y).into();
    let denom = config.atol + config.rtol * y_norm_sq.sqrt();
    let err_scalar = err_norm_sq.sqrt() / denom;

    // Step size update (PI controller simplified to I-only)
    let h_factor = if err_scalar == 0.0 {
        10.0
    } else {
        f64::clamp(0.9 * err_scalar.powf(-0.2), 0.2, 10.0)
    };
    let h_new = dt * h_factor;

    if err_scalar < 1.0 {
        StepResult::Accept { h_new }
    } else {
        StepResult::Reject { h_new }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::prelude::{config::fxhash::ByteF64, PauliSum};
    use crate::lindblad::{LindbladOp, RateMatrix, rhs, rhs_into};
    use crate::solve::SolverCache;

    fn sum1(terms: &[(&str, f64)]) -> PauliSum<ByteF64<1>> {
        let mut s: PauliSum<ByteF64<1>> = PauliSum::builder().n_qubits(1).build();
        for &(w, c) in terms {
            s += (w, c);
        }
        s
    }

    fn get_coeff(s: &PauliSum<ByteF64<1>>, word: &str) -> f64 {
        use ppvm_runtime::prelude::{Trace, PauliWord};
        let w = PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(word);
        s.data().trace(&w)
    }

    fn empty_lindblad() -> LindbladOp<ByteF64<1>> {
        LindbladOp::new(vec![], RateMatrix::from(vec![]))
    }

    #[test]
    fn step_zero_rhs_always_accepted() {
        // With no Hamiltonian and empty Lindblad, rhs = 0.
        // k1 = 0, all k_i = 0, y_new = y, error = 0 => always Accept.
        let y = sum1(&[("X", 1.0)]);
        let lindblad = empty_lindblad();
        let mut cache = SolverCache::new(&y);
        rhs_into(None, &lindblad, &y, &mut cache.k[0]);
        let config = SolverConfig::default();
        match step(None, &lindblad, &y, 0.1, &config, &mut cache) {
            StepResult::Accept { .. } => {
                // y_new is in cache.y_scratch after Accept
                assert!((get_coeff(&cache.y_scratch, "X") - 1.0).abs() < 1e-15);
                assert_eq!(get_coeff(&cache.y_scratch, "Y"), 0.0);
            }
            StepResult::Reject { .. } => panic!("expected Accept"),
        }
    }

    #[test]
    fn step_larmor_small_dt() {
        // H = 0.5*Z, P(0) = X, dt = 0.01
        // Exact: P(t) = cos(t)*X - sin(t)*Y
        // First-order Taylor: X - 0.01*Y, error O(dt^2) ≈ 1e-4
        let h = sum1(&[("Z", 0.5)]);
        let y = sum1(&[("X", 1.0)]);
        let lindblad = empty_lindblad();
        let mut cache = SolverCache::new(&y);
        rhs_into(Some(&h), &lindblad, &y, &mut cache.k[0]);
        let config = SolverConfig::default();
        match step(Some(&h), &lindblad, &y, 0.01, &config, &mut cache) {
            StepResult::Accept { .. } => {
                // y_new is in cache.y_scratch
                let x_coeff = get_coeff(&cache.y_scratch, "X");
                let y_coeff = get_coeff(&cache.y_scratch, "Y");
                assert!((x_coeff - 1.0).abs() < 1e-4, "X coeff: {x_coeff}");
                assert!((y_coeff - (-0.01)).abs() < 1e-4, "Y coeff: {y_coeff}");
            }
            StepResult::Reject { .. } => panic!("expected Accept"),
        }
    }

    #[test]
    fn step_large_dt_rejected() {
        // H = 0.5*Z, P = X, very large dt with default tolerances.
        // The local error grows as O(dt^5), so a large step must be rejected.
        let h = sum1(&[("Z", 0.5)]);
        let y = sum1(&[("X", 1.0)]);
        let lindblad = empty_lindblad();
        let mut cache = SolverCache::new(&y);
        rhs_into(Some(&h), &lindblad, &y, &mut cache.k[0]);
        // Very tight tolerances to ensure rejection even for moderate dt
        let config = SolverConfig { rtol: 1e-10, atol: 1e-12, ..SolverConfig::default() };
        match step(Some(&h), &lindblad, &y, 5.0, &config, &mut cache) {
            StepResult::Reject { h_new } => {
                assert!(h_new < 5.0, "h_new should be smaller: {h_new}");
            }
            StepResult::Accept { .. } => panic!("expected Reject"),
        }
    }

    #[test]
    fn step_fsal_k_next_equals_rhs_y_new() {
        // After Accept, cache.k[6] must equal rhs(cache.y_scratch) computed independently.
        let h = sum1(&[("Z", 0.5)]);
        let y = sum1(&[("X", 1.0)]);
        let lindblad = empty_lindblad();
        let mut cache = SolverCache::new(&y);
        rhs_into(Some(&h), &lindblad, &y, &mut cache.k[0]);
        let config = SolverConfig::default();
        match step(Some(&h), &lindblad, &y, 0.01, &config, &mut cache) {
            StepResult::Accept { .. } => {
                let k7_independent = rhs(Some(&h), &lindblad, &cache.y_scratch);
                for w in &["X", "Y", "Z", "I"] {
                    let a = get_coeff(&cache.k[6], w);
                    let b = get_coeff(&k7_independent, w);
                    assert!((a - b).abs() < 1e-15, "k_next differs at {w}: {a} vs {b}");
                }
            }
            StepResult::Reject { .. } => panic!("expected Accept"),
        }
    }

    // ---- Task 8 tests ----

    // ---- Task 15 tests ----

    #[test]
    fn truncate_state_does_not_accumulate() {
        // Verify that y_scratch is truncated after each accepted step so sub-threshold
        // Pauli strings do not accumulate. We drive a single step directly so the test
        // does not depend on the adaptive step-size controller.
        //
        // Setup: 1-qubit X-dephasing, P(0) = Z, CoefficientThreshold(0.5).
        // After one step the only surviving term should be Z (with coefficient ~1).
        // Without truncation the step would also insert tiny I, X, Y terms from the
        // k-vector contributions, growing the state beyond 1 entry.
        use ppvm_runtime::strategy::CoefficientThreshold;
        use ppvm_runtime::config::fxhash::ByteF64;
        use crate::lindblad::{CollapseOp, JumpOp, LindbladOp, RateMatrix};
        use crate::solve::SolverCache;
        use ppvm_runtime::prelude::{PauliWord, PhasedPauliWord, PauliSum};

        type S = ByteF64<1, CoefficientThreshold>;

        let ppw_s = |pauli: &str, phase: u8|
            -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher,
                               PauliWord<[u8; 1], fxhash::FxBuildHasher>>
        {
            PhasedPauliWord::build_from_word(
                PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli), phase)
        };

        let mut c = CollapseOp::<S>::new(1);
        c.push(ppw_s("X", 0), 1.0);  // c = X (dephasing)
        let lindblad = LindbladOp::new(vec![JumpOp::Generic(c)], RateMatrix::from(vec![1.0]));

        // Use an aggressive threshold so anything below 0.5 is truncated.
        let strat = CoefficientThreshold(0.5);
        let mut y: PauliSum<S> = PauliSum::builder().n_qubits(1).strategy(strat).build();
        y += ("Z", 1.0_f64);

        let config = crate::solve::SolverConfig {
            h0: Some(0.01),  // fixed step so the test is deterministic
            ..crate::solve::SolverConfig::default()
        };
        let mut cache = SolverCache::new(&y);
        rhs_into(None, &lindblad, &y, &mut cache.k[0]);

        match step(None, &lindblad, &y, 0.01, &config, &mut cache) {
            StepResult::Accept { .. } => {
                // y_scratch now holds the truncated new state.
                // Z decays but remains >> 0.5; I/X/Y contributions are O(dt) << 0.5
                // and must have been truncated away.
                let size = cache.y_scratch.data().len();
                assert_eq!(
                    size, 1,
                    "expected only Z to survive truncation, got {size} entries"
                );
            }
            StepResult::Reject { .. } => panic!("expected Accept for dt=0.01"),
        }
    }

    #[test]
    fn truncate_state_preserves_accuracy() {
        // Spontaneous emission: c = X + iY, γ = 1, P(0) = Z, no Hamiltonian.
        // Analytic: <Z>(t) = exp(-8t).  Verify the result still matches after adding
        // the truncation step (regression guard).
        use ppvm_runtime::prelude::{PauliWord, PhasedPauliWord, Trace};
        use ppvm_runtime::strategy::CoefficientThreshold;
        use ppvm_runtime::config::fxhash::ByteF64;
        use crate::lindblad::{CollapseOp, JumpOp, LindbladOp, RateMatrix};
        use crate::solve::solve;

        type S = ByteF64<1, CoefficientThreshold>;

        let ppw = |pauli: &str, phase: u8|
            -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher,
                               PauliWord<[u8; 1], fxhash::FxBuildHasher>>
        {
            PhasedPauliWord::build_from_word(
                PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli), phase)
        };

        let mut c = CollapseOp::<S>::new(1);
        c.push(ppw("X", 0), 1.0);
        c.push(ppw("Y", 1), 1.0);
        let lindblad = LindbladOp::new(vec![JumpOp::Generic(c)], RateMatrix::from(vec![1.0]));

        let strat = CoefficientThreshold(1e-6);
        let mut initial: ppvm_runtime::prelude::PauliSum<S> =
            ppvm_runtime::prelude::PauliSum::builder().n_qubits(1).strategy(strat).build();
        initial += ("Z", 1.0_f64);

        let save_at = [0.25, 0.5, 1.0, 2.0];
        let config = crate::solve::SolverConfig::default();
        let (ts, rs) = solve(
            None, &lindblad, &initial,
            (0.0, 2.0), &save_at,
            |_, p| {
                let w = PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from("Z");
                p.data().trace(&w)
            },
            config,
        );

        assert_eq!(ts.as_slice(), &save_at);
        for (t, r) in ts.iter().zip(rs.iter()) {
            let expected = (-8.0 * t).exp();
            assert!(
                (r - expected).abs() < 1e-4,
                "at t={t}: got {r}, expected {expected}"
            );
        }
    }

    #[test]
    fn estimate_h0_uses_specified_h0() {
        // If config.h0 = Some(x), returns x regardless of the system.
        let y = sum1(&[("X", 1.0)]);
        let h = sum1(&[("Z", 0.5)]);
        let config = SolverConfig { h0: Some(0.1), ..SolverConfig::default() };
        let h = estimate_h0(Some(&h), &empty_lindblad(), &y, (0.0, 1.0), &config);
        assert!((h - 0.1).abs() < 1e-15);
    }

    #[test]
    fn estimate_h0_zero_rhs_fallback() {
        // With ham=None and empty Lindblad, f0 = 0 so d1 = 0 => fallback h0 = 1e-6.
        let y = sum1(&[("X", 1.0)]);
        let config = SolverConfig::default();
        let h = estimate_h0(None, &empty_lindblad(), &y, (0.0, 1.0), &config);
        // d1 = 0 < 1e-5, so h0 = 1e-6; h1 also falls back; result ≥ hmin
        assert!(h > 0.0);
        assert!(h <= 1.0); // within t_span
    }

    #[test]
    fn estimate_h0_nontrivial_system() {
        // H = 0.5*Z, P = X. Estimated h0 is positive, finite, ≤ t_span length.
        let h = sum1(&[("Z", 0.5)]);
        let y = sum1(&[("X", 1.0)]);
        let config = SolverConfig::default();
        let h0 = estimate_h0(Some(&h), &empty_lindblad(), &y, (0.0, 1.0), &config);
        assert!(h0 > 0.0, "h0 should be positive");
        assert!(h0.is_finite(), "h0 should be finite");
        assert!(h0 <= 1.0, "h0 should be ≤ t_span length");
    }
}
