/// Superradiance example — Budget truncation showcase.
///
/// Demonstrates that a tight `Budget` cap achieves a meaningful speedup over a
/// generous baseline while keeping observable error small.  Per-stage truncation
/// (Task 26) makes both variants work with `SolverConfig::default()` — no manual
/// rtol tuning is required.
///
/// Variants
/// --------
/// 1. Baseline : Budget { target=2000 } — generous cap; state never hits the limit,
///               so this behaves like untruncated evolution.
/// 2. Budget   : Budget { target=200  } — tight cap; DOPRI5 takes fewer RHS evaluations
///               per unit time because the state is smaller.  Accuracy loss is small.
///
/// System: n=5, 25 Lindblad raising terms (< PAR_THRESHOLD=200 → sequential path).
/// Change N to 15 to engage Rayon parallelism (225 terms > threshold).
///
/// Run with:  cargo run --example superradiance --release
use std::time::Instant;

use ppvm_runtime::{
    config::indexmap::ByteFxHashF64,
    prelude::*,
};
use ppvm_timeevolve::{Budget, JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix, SolverConfig, solve::solve};

const N:      usize = 5;
const NBYTES: usize = 1;
const GAMMA0: f64   = 1.0;
const D:      f64   = 0.1;
const TMAX:   f64   = 0.5;
const TSTEPS: usize = 20;

/// Baseline: generous cap — state stays well below this for n=5.
const BASE_TARGET:   usize = 2000;
/// Budget: tight cap — forces truncation once the state grows.
const BUDGET_TARGET: usize = 200;

type SB = ByteFxHashF64<NBYTES, Budget>;

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

fn build_ops() -> LindbladOp<SB> {
    let ops: Vec<JumpOp<SB>> = (0..N)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Raise }))
        .collect();
    LindbladOp::new(ops, rate_matrix())
}

/// Initial observable O(0). Propagated under `dO/dt = L†(O)` — NOT the density matrix.
fn initial_observable(strat: Budget) -> PauliSum<SB> {
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

    // ── Variant 1: Baseline ──────────────────────────────────────────────────
    let t0 = Instant::now();
    let (_, base_out) = solve(
        None, &lindblad, &initial_observable(Budget { target: BASE_TARGET }),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| (p.trace(&pattern), p.data().len()),
        SolverConfig::default(),
    );
    let t_base = t0.elapsed();

    // ── Variant 2: Budget (tight cap, default rtol) ──────────────────────────
    // Per-stage truncation keeps the DOPRI5 error estimate consistent with the
    // truncated ODE, so no rtol adjustment is needed.
    let t1 = Instant::now();
    let (_, bud_out) = solve(
        None, &lindblad, &initial_observable(Budget { target: BUDGET_TARGET }),
        (0.0, TMAX), &save_at,
        |_, p: &PauliSum<SB>| (p.trace(&pattern), p.data().len()),
        SolverConfig::default(),
    );
    let t_bud = t1.elapsed();

    // ── Fidelity (max |variant − baseline| over all save points) ────────────
    let fidelity = |out: &[(f64, usize)]| -> f64 {
        base_out.iter().zip(out.iter())
            .map(|((b, _), (v, _))| (b - v).abs())
            .fold(0.0_f64, f64::max)
    };

    let final_p = |out: &[(f64, usize)]| out.last().map(|&(_, sz)| sz).unwrap_or(0);
    let spd_bud  = t_base.as_secs_f64() / t_bud.as_secs_f64();

    // ── Print table ──────────────────────────────────────────────────────────
    println!("n={N}  tmax={TMAX}  steps={TSTEPS}  Lindblad terms={}", N * N);
    println!("(sequential path: {} terms < PAR_THRESHOLD=200; change N to 15 to engage Rayon)", N * N);
    println!();
    println!("{:<20} {:>10} {:>10} {:>10} {:>18}",
        "Variant", "Wall time", "Speedup", "|P| final", "Max fidelity err");
    println!("{}", "─".repeat(72));
    println!("{:<20} {:>10.3?} {:>10} {:>10} {:>18}",
        "Baseline", t_base, "1.00×", final_p(&base_out), "—");
    println!("{:<20} {:>10.3?} {:>9.2}× {:>10} {:>18.2e}",
        "Budget", t_bud, spd_bud, final_p(&bud_out), fidelity(&bud_out));
    println!();
    println!("Baseline : Budget{{target={BASE_TARGET}}}, default rtol — generous cap, untruncated behaviour");
    println!("Budget   : Budget{{target={BUDGET_TARGET}}}, default rtol — tight cap, per-stage truncation keeps DOPRI5 stable");
}
