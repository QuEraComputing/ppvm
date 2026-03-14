use std::ops::{Mul, MulAssign};

use ppvm_runtime::prelude::{
    ACMapAddAssign, ACMapIter, Config, PauliSum, PhasedPauliWord, Trace,
};

use crate::lindblad::{LindbladOp, rhs};
use crate::solve::SolverConfig;

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

#[allow(dead_code)] // fields used from solve.rs (Task 9)
pub(crate) enum StepResult<T: Config> {
    Accept { y_new: PauliSum<T>, k_next: PauliSum<T>, h_new: f64 },
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
        *target += (w.clone(), s * *c);
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
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign
        + Copy
        + std::ops::Mul<Output = T::Coeff>
        + std::iter::Sum
        + Into<f64>,
    T::PauliWordType: Clone,
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
/// `k1` is `rhs(y)` from the previous step (FSAL).
/// Returns `Accept` if the local error is within tolerance, `Reject` otherwise.
/// In either case `h_new` gives the suggested next step size.
#[allow(dead_code)] // called from solve.rs (Task 9)
pub(crate) fn step<T: Config>(
    ham: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    y: &PauliSum<T>,
    k1: PauliSum<T>,
    dt: f64,
    config: &SolverConfig,
) -> StepResult<T>
where
    for<'a> T::Map: ACMapIter<'a, Item = (&'a T::PauliWordType, &'a T::Coeff)>,
    T::Map: ACMapAddAssign<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType>,
    T::Coeff: std::ops::AddAssign
        + Copy
        + std::ops::Mul<Output = T::Coeff>
        + std::iter::Sum
        + Into<f64>,
    T::PauliWordType: Clone,
    PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        Mul<Output = PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + MulAssign
        + Clone,
    f64: Into<T::Coeff>,
    for<'a> T::Map: Trace<'a, T::PauliWordType, Output = T::Coeff>,
{
    // Stage 2
    let k2 = {
        let mut y2 = y.clone();
        add_scaled(&mut y2, &k1, dt * A21);
        rhs(ham, lindblad, &y2)
    };

    // Stage 3
    let k3 = {
        let mut y3 = y.clone();
        add_scaled(&mut y3, &k1, dt * A31);
        add_scaled(&mut y3, &k2, dt * A32);
        rhs(ham, lindblad, &y3)
    };

    // Stage 4
    let k4 = {
        let mut y4 = y.clone();
        add_scaled(&mut y4, &k1, dt * A41);
        add_scaled(&mut y4, &k2, dt * A42);
        add_scaled(&mut y4, &k3, dt * A43);
        rhs(ham, lindblad, &y4)
    };

    // Stage 5
    let k5 = {
        let mut y5 = y.clone();
        add_scaled(&mut y5, &k1, dt * A51);
        add_scaled(&mut y5, &k2, dt * A52);
        add_scaled(&mut y5, &k3, dt * A53);
        add_scaled(&mut y5, &k4, dt * A54);
        rhs(ham, lindblad, &y5)
    };

    // Stage 6
    let k6 = {
        let mut y6 = y.clone();
        add_scaled(&mut y6, &k1, dt * A61);
        add_scaled(&mut y6, &k2, dt * A62);
        add_scaled(&mut y6, &k3, dt * A63);
        add_scaled(&mut y6, &k4, dt * A64);
        add_scaled(&mut y6, &k5, dt * A65);
        rhs(ham, lindblad, &y6)
    };

    // 5th-order update: y_new = y + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
    // b_2 = 0 and b_7 = 0 (FSAL: k7 not used in y_new)
    let mut y_new = y.clone();
    add_scaled(&mut y_new, &k1, dt * B1);
    add_scaled(&mut y_new, &k3, dt * B3);
    add_scaled(&mut y_new, &k4, dt * B4);
    add_scaled(&mut y_new, &k5, dt * B5);
    add_scaled(&mut y_new, &k6, dt * B6);

    // Stage 7 (FSAL: k_next = k1 of the next step, also used in error estimate)
    let k7 = rhs(ham, lindblad, &y_new);

    // Error estimate: e = dt * (e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
    // e_2 = 0
    let mut err_vec = PauliSum::<T>::builder().n_qubits(y.n_qubits()).build();
    add_scaled(&mut err_vec, &k1, dt * E1);
    add_scaled(&mut err_vec, &k3, dt * E3);
    add_scaled(&mut err_vec, &k4, dt * E4);
    add_scaled(&mut err_vec, &k5, dt * E5);
    add_scaled(&mut err_vec, &k6, dt * E6);
    add_scaled(&mut err_vec, &k7, dt * E7);

    // Error norm: err = ||e|| / (atol + rtol * ||y||)
    let err_norm_sq: f64 = err_vec.overlap(&err_vec).into();
    let y_norm_sq: f64 = y.overlap(y).into();
    let denom = config.atol + config.rtol * y_norm_sq.sqrt();
    let err = err_norm_sq.sqrt() / denom;

    // Step size update (PI controller simplified to I-only)
    let h_factor = if err == 0.0 {
        10.0
    } else {
        f64::clamp(0.9 * err.powf(-0.2), 0.2, 10.0)
    };
    let h_new = dt * h_factor;

    if err < 1.0 {
        StepResult::Accept { y_new, k_next: k7, h_new }
    } else {
        StepResult::Reject { h_new }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::prelude::{config::fxhash::ByteF64, PauliSum};
    use crate::lindblad::{LindbladOp, RateMatrix};

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
        let k1 = rhs(None, &empty_lindblad(), &y);
        let config = SolverConfig::default();
        match step(None, &empty_lindblad(), &y, k1, 0.1, &config) {
            StepResult::Accept { y_new, .. } => {
                assert!((get_coeff(&y_new, "X") - 1.0).abs() < 1e-15);
                assert_eq!(get_coeff(&y_new, "Y"), 0.0);
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
        let k1 = rhs(Some(&h), &empty_lindblad(), &y);
        let config = SolverConfig::default();
        match step(Some(&h), &empty_lindblad(), &y, k1, 0.01, &config) {
            StepResult::Accept { y_new, .. } => {
                let x_coeff = get_coeff(&y_new, "X");
                let y_coeff = get_coeff(&y_new, "Y");
                // Accept if within 1e-4 of first-order Taylor
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
        let k1 = rhs(Some(&h), &empty_lindblad(), &y);
        // Very tight tolerances to ensure rejection even for moderate dt
        let config = SolverConfig { rtol: 1e-10, atol: 1e-12, ..SolverConfig::default() };
        match step(Some(&h), &empty_lindblad(), &y, k1, 5.0, &config) {
            StepResult::Reject { h_new } => {
                assert!(h_new < 5.0, "h_new should be smaller: {h_new}");
            }
            StepResult::Accept { .. } => panic!("expected Reject"),
        }
    }

    #[test]
    fn step_fsal_k_next_equals_rhs_y_new() {
        // k_next from Accept must equal rhs(y_new) computed independently.
        let h = sum1(&[("Z", 0.5)]);
        let y = sum1(&[("X", 1.0)]);
        let k1 = rhs(Some(&h), &empty_lindblad(), &y);
        let config = SolverConfig::default();
        let lindblad = empty_lindblad();
        match step(Some(&h), &lindblad, &y, k1, 0.01, &config) {
            StepResult::Accept { y_new, k_next, .. } => {
                let k7_independent = rhs(Some(&h), &lindblad, &y_new);
                // k_next should equal rhs(y_new) to floating-point precision
                for w in &["X", "Y", "Z", "I"] {
                    let a = get_coeff(&k_next, w);
                    let b = get_coeff(&k7_independent, w);
                    assert!((a - b).abs() < 1e-15, "k_next differs at {w}: {a} vs {b}");
                }
            }
            StepResult::Reject { .. } => panic!("expected Accept"),
        }
    }

    // ---- Task 8 tests ----

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
