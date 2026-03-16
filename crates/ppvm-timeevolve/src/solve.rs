use std::ops::{Mul, MulAssign};

use ppvm_runtime::prelude::{
    ACMapAddAssign, ACMapIter, Config, PauliSum, PhasedPauliWord, Trace,
};

use crate::dopri5::{StepResult, estimate_h0, step};
use crate::lindblad::{LindbladOp, rhs};

pub struct SolverConfig {
    pub rtol:  f64,
    pub atol:  f64,
    pub h0:    Option<f64>,
    pub hmin:  f64,
    pub hmax:  f64,
}

impl Default for SolverConfig {
    fn default() -> Self {
        SolverConfig {
            rtol:  1e-6,
            atol:  1e-9,
            h0:    None,
            hmin:  1e-12,
            hmax:  f64::INFINITY,
        }
    }
}

/// Advance `state` in-place from `t_span.0` to `t_span.1`.
///
/// At each time in `save_at` (sorted, within t_span), the step is capped
/// exactly to that time and `callback` is invoked.  Returns the recorded
/// times (equal to `save_at`) and the corresponding callback outputs.
pub fn solve_mut<T: Config, R, F>(
    hamiltonian: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    state: &mut PauliSum<T>,
    t_span: (f64, f64),
    save_at: &[f64],
    callback: F,
    config: SolverConfig,
) -> (Vec<f64>, Vec<R>)
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
    F: Fn(f64, &PauliSum<T>) -> R,
{
    let mut t = t_span.0;
    let mut dt = estimate_h0(hamiltonian, lindblad, state, t_span, &config);
    let mut k1 = rhs(hamiltonian, lindblad, state);
    let mut times = Vec::with_capacity(save_at.len());
    let mut results = Vec::with_capacity(save_at.len());

    for &t_save in save_at {
        // Advance from t to t_save using adaptive steps.
        loop {
            let remaining = t_save - t;
            // Stop once we are negligibly close to t_save.
            if remaining <= f64::EPSILON * t_save.abs().max(1.0) {
                break;
            }
            let dt_capped = dt.min(remaining);

            match step(hamiltonian, lindblad, state, k1.clone(), dt_capped, &config) {
                StepResult::Accept { y_new, k_next, h_new } => {
                    t += dt_capped;
                    *state = y_new;
                    k1 = k_next;
                    dt = h_new;
                }
                StepResult::Reject { h_new } => {
                    dt = h_new.max(config.hmin);
                }
            }
        }

        times.push(t_save);
        results.push(callback(t_save, state));
    }

    (times, results)
}

/// Clone `initial` and delegate to [`solve_mut`].  `initial` is not modified.
pub fn solve<T: Config, R, F>(
    hamiltonian: Option<&PauliSum<T>>,
    lindblad: &LindbladOp<T>,
    initial: &PauliSum<T>,
    t_span: (f64, f64),
    save_at: &[f64],
    callback: F,
    config: SolverConfig,
) -> (Vec<f64>, Vec<R>)
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
    F: Fn(f64, &PauliSum<T>) -> R,
{
    let mut state = initial.clone();
    solve_mut(hamiltonian, lindblad, &mut state, t_span, save_at, callback, config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::prelude::{PauliSum, PauliWord, PhasedPauliWord, Trace, config::fxhash::ByteF64};
    use crate::lindblad::{CollapseOp, LindbladOp, RateMatrix};

    type S = ByteF64<1>;

    fn sum1(terms: &[(&str, f64)]) -> PauliSum<S> {
        let mut s: PauliSum<S> = PauliSum::builder().n_qubits(1).build();
        for &(w, c) in terms {
            s += (w, c);
        }
        s
    }

    fn sum2(terms: &[(&str, f64)]) -> PauliSum<S> {
        let mut s: PauliSum<S> = PauliSum::builder().n_qubits(2).build();
        for &(w, c) in terms {
            s += (w, c);
        }
        s
    }

    fn get_coeff(s: &PauliSum<S>, word: &str) -> f64 {
        let w = PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(word);
        s.data().trace(&w)
    }

    fn empty_lindblad() -> LindbladOp<S> {
        LindbladOp::new(vec![], RateMatrix::from(vec![]))
    }

    fn x_sum_overlap(p: &PauliSum<S>) -> f64 {
        // <X|P> = coefficient of X in P
        get_coeff(p, "X")
    }

    #[test]
    fn solve_config_defaults() {
        let c = SolverConfig::default();
        assert_eq!(c.rtol, 1e-6);
        assert_eq!(c.atol, 1e-9);
        assert_eq!(c.h0, None);
        assert_eq!(c.hmin, 1e-12);
        assert_eq!(c.hmax, f64::INFINITY);
    }

    #[test]
    fn solve_empty_save_at() {
        // Empty save_at returns ([], []) without panic.
        let h = sum1(&[("Z", 0.5)]);
        let mut state = sum1(&[("X", 1.0)]);
        let config = SolverConfig::default();
        let (ts, rs) = solve_mut(
            Some(&h), &empty_lindblad(), &mut state,
            (0.0, 1.0), &[], |_, _| 0.0_f64, config,
        );
        assert!(ts.is_empty());
        assert!(rs.is_empty());
    }

    #[test]
    fn solve_larmor_single_save() {
        // H = 0.5*Z, P(0) = X, save at t = 1.0.
        // Exact solution: P(t) = cos(t)*X - sin(t)*Y
        // <X|P(1)> = cos(1) ≈ 0.5403023
        let h = sum1(&[("Z", 0.5)]);
        let initial = sum1(&[("X", 1.0)]);
        let config = SolverConfig::default();
        let (ts, rs) = solve(
            Some(&h), &empty_lindblad(), &initial,
            (0.0, 1.0), &[1.0],
            |_, p| x_sum_overlap(p),
            config,
        );
        assert_eq!(ts, vec![1.0]);
        let expected = (1.0_f64).cos();
        assert!(
            (rs[0] - expected).abs() < 1e-4,
            "got {}, expected {}",
            rs[0], expected
        );
    }

    #[test]
    fn solve_larmor_multiple_saves() {
        // Same as above, but save at 0.25, 0.5, 0.75, 1.0.
        let h = sum1(&[("Z", 0.5)]);
        let initial = sum1(&[("X", 1.0)]);
        let config = SolverConfig::default();
        let save_at = [0.25, 0.5, 0.75, 1.0];
        let (ts, rs) = solve(
            Some(&h), &empty_lindblad(), &initial,
            (0.0, 1.0), &save_at,
            |_, p| x_sum_overlap(p),
            config,
        );
        assert_eq!(ts.as_slice(), &save_at);
        for (t_s, r) in ts.iter().zip(rs.iter()) {
            let expected = t_s.cos();
            assert!(
                (r - expected).abs() < 1e-4,
                "at t={t_s}: got {r}, expected {expected}"
            );
        }
    }

    #[test]
    fn solve_vs_solve_mut_identical() {
        // solve and solve_mut must produce identical results.
        // solve must not modify initial.
        let h = sum1(&[("Z", 0.5)]);
        let initial = sum1(&[("X", 1.0)]);
        let mut state = initial.clone();
        let save_at = [0.5, 1.0];
        let config1 = SolverConfig::default();
        let config2 = SolverConfig::default();

        let (ts1, rs1) = solve_mut(
            Some(&h), &empty_lindblad(), &mut state,
            (0.0, 1.0), &save_at, |_, p| x_sum_overlap(p), config1,
        );
        let (ts2, rs2) = solve(
            Some(&h), &empty_lindblad(), &initial,
            (0.0, 1.0), &save_at, |_, p| x_sum_overlap(p), config2,
        );

        assert_eq!(ts1, ts2);
        for (a, b) in rs1.iter().zip(rs2.iter()) {
            assert!((a - b).abs() < 1e-14, "solve vs solve_mut differ: {a} vs {b}");
        }

        // initial unchanged
        assert!((get_coeff(&initial, "X") - 1.0).abs() < 1e-15);
        assert_eq!(get_coeff(&initial, "Y"), 0.0);
    }

    #[test]
    fn solve_two_qubit_correlated_dephasing() {
        // Two qubits with correlated Z-dephasing (dipole-dipole-type off-diagonal rate).
        // c_1 = ZI (Z on qubit 1, I on qubit 2), c_2 = IZ (I on qubit 1, Z on qubit 2).
        // Rate matrix: Γ = [[γ, γ_12], [γ_12, γ]] with γ=1.0, γ_12=0.5 (γ_12 < γ).
        // No Hamiltonian.  P(0) = XX + YY, 2-qubit PauliSum.
        //
        // Analytic derivation:
        //   c_1†c_1 = ZI·ZI = II,  c_2†c_2 = IZ·IZ = II,  c_1†c_2 = c_2†c_1 = ZZ.
        //
        //   Diagonal term (i=j=1, rate γ): 2 ZI P ZI − {II, P} = 2(ZXZ)⊗X − 2XX = −4XX.
        //   Diagonal term (i=j=2, rate γ): 2 X⊗(ZXZ) − 2XX = −4XX.
        //   → diagonal total on XX: −8γ XX.
        //
        //   Off-diagonal (1,2) and (2,1), rate γ_12 each:
        //     sandwich ZI (XX) IZ = (ZX)⊗(XI)(IZ) → (iY)⊗(−iY) = YY  →  +2γ_12 YY each.
        //     anticommutator {ZZ, XX}: (ZX⊗ZX)+(XZ⊗XZ) = (iY)(iY)+(−iY)(−iY) = −2YY
        //     → each off-diagonal term contributes 2γ_12(YY−(−YY)) = 4γ_12 YY.
        //     Two pairs total: 8γ_12 YY.
        //   → L†(XX) = −8γ XX + 8γ_12 YY.
        //   → L†(YY) = −8γ YY + 8γ_12 XX.   (by X↔Y symmetry)
        //
        //   Normal modes XX±YY:
        //     d(XX+YY)/dt = −8(γ−γ_12)(XX+YY)   eigenvalue λ₊ = −8(γ−γ_12)
        //     d(XX−YY)/dt = −8(γ+γ_12)(XX−YY)   eigenvalue λ₋ = −8(γ+γ_12)
        //
        //   With γ=1, γ_12=0.5: λ₊ = −4,  λ₋ = −12.
        //   P(0) = XX+YY lies purely in the λ₊ mode:
        //     coefficient of XX at time t = e^{−4t}.

        let ppw = |pauli: &str, phase: u8|
            -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, PauliWord<[u8; 1], fxhash::FxBuildHasher>>
        {
            PhasedPauliWord::build_from_word(
                PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli), phase)
        };

        let mut c1 = CollapseOp::<S>::new(2);
        c1.push(ppw("ZI", 0), 1.0);
        let mut c2 = CollapseOp::<S>::new(2);
        c2.push(ppw("IZ", 0), 1.0);
        let lindblad = LindbladOp::new(
            vec![c1, c2],
            RateMatrix::Dense(vec![vec![1.0, 0.5], vec![0.5, 1.0]]),
        );

        let initial = sum2(&[("XX", 1.0), ("YY", 1.0)]);
        let save_at = [0.1, 0.25, 0.5, 1.0];
        let config = SolverConfig::default();
        let (ts, rs) = solve(
            None, &lindblad, &initial,
            (0.0, 1.0), &save_at,
            |_, p| get_coeff(p, "XX"),
            config,
        );

        assert_eq!(ts.as_slice(), &save_at);
        for (t_s, r) in ts.iter().zip(rs.iter()) {
            let expected = (-4.0 * t_s).exp();
            assert!(
                (r - expected).abs() < 1e-4,
                "at t={t_s}: got {r}, expected {expected}"
            );
        }
    }

    #[test]
    fn solve_spontaneous_emission() {
        // Setup: c = X + iY (un-normalised lowering operator), γ = 1, no Hamiltonian.
        // P(0) = Z, single qubit.
        //
        // Analytic derivation:
        //   From Task 5: L(Z) = 8I − 8Z  for c = X+iY, γ = 1.
        //   L(I) = 0 (identity is preserved by the adjoint Lindblad).
        //   Write P(t) = a_I(t)·I + a_Z(t)·Z.  With P(0) = Z: a_I(0)=0, a_Z(0)=1.
        //   dP/dt = L(P)  gives:
        //     d(a_Z)/dt = −8·a_Z   =>   a_Z(t) = e^{−8t}
        //     d(a_I)/dt =  8·a_Z   =>   a_I(t) = 1 − e^{−8t}
        //   Callback extracts coefficient of Z, so expected value = e^{−8t}.

        // Build c = X + iY: X has phase 0, iY is Y with phase 1 (i^1 = i).
        let ppw = |pauli: &str, phase: u8| -> PhasedPauliWord<[u8; 1], fxhash::FxBuildHasher, PauliWord<[u8; 1], fxhash::FxBuildHasher>> {
            PhasedPauliWord::build_from_word(PauliWord::<[u8; 1], fxhash::FxBuildHasher>::from(pauli), phase)
        };
        let mut c = CollapseOp::<S>::new(1);
        c.push(ppw("X", 0), 1.0);
        c.push(ppw("Y", 1), 1.0);
        let lindblad = LindbladOp::new(vec![c], RateMatrix::from(vec![1.0]));

        let initial = sum1(&[("Z", 1.0)]);
        let save_at = [0.25, 0.5, 1.0, 2.0];
        let config = SolverConfig::default();
        let (ts, rs) = solve(
            None, &lindblad, &initial,
            (0.0, 2.0), &save_at,
            |_, p| get_coeff(p, "Z"),
            config,
        );

        assert_eq!(ts.as_slice(), &save_at);
        for (t_s, r) in ts.iter().zip(rs.iter()) {
            let expected = (-8.0 * t_s).exp();
            assert!(
                (r - expected).abs() < 1e-4,
                "at t={t_s}: got {r}, expected {expected}"
            );
        }
    }
}
