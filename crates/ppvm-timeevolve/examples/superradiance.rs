/// Superradiance example — three-variant comparison (Budget truncation showcase).
///
/// Demonstrates that the speedup from `Budget` comes from the *rtol coupling*,
/// not from Rayon.  Rayon parallelism is a separate, additive benefit that
/// engages automatically for n≥8 (256 terms > PAR_THRESHOLD=200).
///
/// Variants
/// --------
/// 1. Baseline      : Budget{target=2000, min_threshold=1e-6}, default rtol=1e-6
///                    — generous cap; acts like CoefficientThreshold(1e-6) here
/// 2. Budget        : Budget{target=200, min_threshold=1e-4}, default rtol=1e-6
///                    — tight cap + looser threshold, but rtol is mismatched:
///                      DOPRI5 fights the k[6] truncation perturbation → slower
/// 3. Budget+rtol   : same Budget, but rtol=10·min_threshold=1e-3
///                    — rtol absorbs the k[6] perturbation (≈ h·|e7|·||L||·||Δy||),
///                      letting DOPRI5 take large steps.  Speedup is visible even
///                      on the sequential path (n=5); at n≥8 Rayon adds more.
///
/// System: n=5, 100 Lindblad terms (< PAR_THRESHOLD=200 → sequential path).
/// Rayon parallelism engages automatically for n≥8 (256 terms > threshold).
///
/// Run with:  cargo run --example superradiance --release
use std::time::Instant;

use ppvm_runtime::{
    config::fxhash::ByteF64,
    prelude::*,
};
use ppvm_timeevolve::{Budget, CollapseOp, JumpOp, LindbladOp, RateMatrix, SolverConfig, solve::solve};

const N: usize = 5;
const NBYTES: usize = 1;
const GAMMA0: f64 = 1.0;
const D: f64 = 0.1;
const TMAX: f64 = 0.5;
const TSTEPS: usize = 20;
/// Baseline: generous cap — acts like CoefficientThreshold(1e-6) for |P| < 2000.
const BASE_TARGET:     usize = 2000;
const BASE_THRESHOLD:  f64   = 1e-6;
/// Budget variants: tight cap + looser threshold.
const BUDGET_TARGET:    usize = 200;
const BUDGET_THRESHOLD: f64   = 1e-4;
/// Practical rtol for Budget+rtol: larger than min_threshold to account for the
/// truncation-induced perturbation of k[6] (≈ h·|e7|·||L||·||Δy||).
/// Rule of thumb: rtol ≈ 10 · min_threshold absorbs that noise and lets DOPRI5
/// take large steps; ODE accuracy then matches observable-level truncation error.
const BUDGET_RTOL: f64 = BUDGET_THRESHOLD * 10.0; // = 1e-3

type W  = PauliWord<[u8; NBYTES], fxhash::FxBuildHasher>;
type SB = ByteF64<NBYTES, Budget>;

fn rate_matrix() -> RateMatrix {
    RateMatrix::Dense(
        (0..N)
            .map(|i| {
                (0..N)
                    .map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs()))
                    .collect()
            })
            .collect(),
    )
}

fn ppw(s: &str, phase: u8) -> PhasedPauliWord<[u8; NBYTES], fxhash::FxBuildHasher, W> {
    PhasedPauliWord::build_from_word(W::from(s), phase)
}

fn build_ops() -> LindbladOp<SB> {
    let t = vec!['I'; N];
    let mut ops: Vec<CollapseOp<SB>> = Vec::with_capacity(N);
    for i in 0..N {
        let mut op = CollapseOp::new(N);
        let mut px = t.clone(); px[i] = 'X';
        let mut py = t.clone(); py[i] = 'Y';
        op.push(ppw(&px.iter().collect::<String>(), 0), 1.0);
        op.push(ppw(&py.iter().collect::<String>(), 3), 1.0);
        ops.push(op);
    }
    LindbladOp::new(ops.into_iter().map(JumpOp::Generic).collect(), rate_matrix())
}

fn initial_state(strat: Budget) -> PauliSum<SB> {
    let mut p = PauliSum::builder().n_qubits(N).strategy(strat).build();
    let t = vec!['I'; N];
    for i in 0..N {
        let mut zi = t.clone(); zi[i] = 'Z';
        p += (zi.iter().collect::<String>(), 1.0_f64);
    }
    p
}

fn main() {
    let lindblad = build_ops();
    let save_at: Vec<f64> = (1..=TSTEPS).map(|i| i as f64 * TMAX / TSTEPS as f64).collect();
    let pattern: PauliPattern = "Z?*".into();

    // Each solve callback returns (trace_value, state_size) in one pass.

    // ── Variant 1: Baseline ──────────────────────────────────────────────────
    let strat_base = Budget { target: BASE_TARGET, min_threshold: BASE_THRESHOLD };
    let t0 = Instant::now();
    let (_, base_out) = solve(
        None, &lindblad, &initial_state(strat_base),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| (p.trace(&pattern), p.data().len()),
        SolverConfig::default(),
    );
    let t_base = t0.elapsed();

    // ── Variant 2: Budget (default rtol — mismatched) ────────────────────────
    let strat_bud = Budget { target: BUDGET_TARGET, min_threshold: BUDGET_THRESHOLD };
    let t1 = Instant::now();
    let (_, bud_out) = solve(
        None, &lindblad, &initial_state(strat_bud),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| (p.trace(&pattern), p.data().len()),
        SolverConfig::default(),
    );
    let t_bud = t1.elapsed();

    // ── Variant 3: Budget+rtol (matched rtol — sequential) ──────────────────
    // The speedup here comes entirely from the rtol adjustment, not from Rayon.
    // At n≥8 (256 terms > PAR_THRESHOLD=200), Rayon would add further speedup.
    let t2 = Instant::now();
    let (_, bud_rtol_out) = solve(
        None, &lindblad, &initial_state(strat_bud),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| (p.trace(&pattern), p.data().len()),
        // rtol ≈ 10 · min_threshold absorbs the k[6] perturbation from Budget
        // truncation, preventing step rejections while still allowing large steps.
        SolverConfig { rtol: BUDGET_RTOL, ..SolverConfig::default() },
    );
    let t_bud_rtol = t2.elapsed();

    // ── Fidelity (max |variant − baseline| over all save points) ────────────
    let fidelity = |out: &[(f64, usize)]| -> f64 {
        base_out.iter().zip(out.iter())
            .map(|((b, _), (v, _))| (b - v).abs())
            .fold(0.0_f64, f64::max)
    };

    let final_p = |out: &[(f64, usize)]| out.last().map(|&(_, sz)| sz).unwrap_or(0);

    let spd_bud      = t_base.as_secs_f64() / t_bud.as_secs_f64();
    let spd_bud_rtol = t_base.as_secs_f64() / t_bud_rtol.as_secs_f64();

    // ── Print table ──────────────────────────────────────────────────────────
    println!("n={N}  tmax={TMAX}  steps={TSTEPS}  Lindblad terms={}", 4 * N * N);
    println!("(sequential path: {} terms < PAR_THRESHOLD=200; change N to 8 to engage Rayon)", 4 * N * N);
    println!();
    println!("{:<20} {:>10} {:>10} {:>10} {:>18}",
        "Variant", "Wall time", "Speedup", "|P| final", "Max fidelity err");
    println!("{}", "─".repeat(72));
    println!("{:<20} {:>10.3?} {:>10} {:>10} {:>18}",
        "Baseline", t_base, "1.00×", final_p(&base_out), "—");
    println!("{:<20} {:>10.3?} {:>9.2}× {:>10} {:>18.2e}",
        "Budget", t_bud, spd_bud, final_p(&bud_out), fidelity(&bud_out));
    println!("{:<20} {:>10.3?} {:>9.2}× {:>10} {:>18.2e}",
        "Budget+rtol", t_bud_rtol, spd_bud_rtol, final_p(&bud_rtol_out), fidelity(&bud_rtol_out));
    println!();
    let default_rtol = SolverConfig::default().rtol;
    println!("{:<20} rtol={:.0e}  Budget{{target={}, min_threshold={:.0e}}}",
        "Baseline:",     default_rtol,    BASE_TARGET,   BASE_THRESHOLD);
    println!("{:<20} rtol={:.0e}  Budget{{target={}, min_threshold={:.0e}}}  ← mismatched",
        "Budget:",        default_rtol,   BUDGET_TARGET, BUDGET_THRESHOLD);
    println!("{:<20} rtol={:.0e}  Budget{{target={}, min_threshold={:.0e}}}  ← matched (10×min_threshold); at n≥8 Rayon also engages",
        "Budget+rtol:",   BUDGET_RTOL,    BUDGET_TARGET, BUDGET_THRESHOLD);
}
