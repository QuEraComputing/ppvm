/// Scaling study: n-atom superradiance with different truncation strategies.
///
/// Investigates the central question: can truncation keep computation polynomial
/// in n while maintaining sufficient accuracy?
///
/// System
/// ------
/// n atoms, each with a lowering operator σ_−.  Rate matrix γ_ij = γ₀/(1+D·|i−j|)
/// (nearest-neighbour-dominated collective decay).  Observable: ⟨Z₀⟩ at T=0.5.
///
/// Strategies
/// ----------
/// Reference        : CoefficientThreshold(1e-8) — treated as exact
/// CT(1e-3)         : CoefficientThreshold(1e-3) — coefficient pruning
/// MPW(2)           : MaxPauliWeight(2)  — O(n²) state size bound
/// MPW(4)           : MaxPauliWeight(4)  — O(n⁴) state size bound
/// Budget(300,1e-4) : hard cap + matched rtol (10×min_threshold) for DOPRI5
///
/// Scaling
/// -------
/// - Reference |P| grows exponentially; wall time grows correspondingly.
/// - MPW(w) caps |P| ≤ Σ_{k=0}^{w} C(n,k)·3^k  (polynomial in n for fixed w).
/// - Budget(target) caps |P| ≤ target regardless of n.
/// - CT(ε) prunes dynamically; |P| depends on how many coefficients exceed ε.
///
/// Run with:  cargo run --example scaling --release
use std::time::Instant;

use ppvm_runtime::{
    config::fxhash::ByteF64,
    prelude::*,
    strategy::{CoefficientThreshold, MaxPauliWeight},
};
use ppvm_timeevolve::{
    Budget, JumpOp, LadderDirection, LadderOp, LindbladOp, RateMatrix, SolverConfig,
    solve::solve,
};

// ── Physical parameters ───────────────────────────────────────────────────────
const GAMMA0: f64 = 1.0;
const D:      f64 = 0.1;
const TMAX:   f64 = 0.5;

/// Atom counts to sweep.  Reference (CoefficientThreshold 1e-8) becomes slow at n≥8
/// because |P| grows exponentially — change to &[2, 4, 6, 8] to see that regime.
const N_VALUES: &[usize] = &[2, 4, 6];

// ── Strategy parameters ───────────────────────────────────────────────────────
const REF_THRESHOLD:   f64   = 1e-8;
const CT_THRESHOLD:    f64   = 1e-3;
/// Matched rtol for CT: same logic as Budget — rtol ≈ 10 × threshold avoids
/// DOPRI5 fighting the coefficient-truncation perturbation on k[6].
const CT_RTOL:         f64   = CT_THRESHOLD * 10.0;
const MPW2_WEIGHT:     usize = 2;
const MPW4_WEIGHT:     usize = 4;
const BUD_TARGET:      usize = 300;
const BUD_MIN_THRESH:  f64   = 1e-4;
/// Matched rtol for Budget: absorbs the k[6] truncation perturbation in DOPRI5.
/// Rule of thumb: rtol ≈ 10 × min_threshold  (see superradiance example for derivation).
const BUD_RTOL:        f64   = BUD_MIN_THRESH * 10.0;

// ── Config type aliases ───────────────────────────────────────────────────────
// 2 bytes → up to 16 qubits (covers N_VALUES max = 8 with room to spare).
type Sct  = ByteF64<2, CoefficientThreshold>;
type Smpw = ByteF64<2, MaxPauliWeight>;
type Sbud = ByteF64<2, Budget>;

// ── Helpers ───────────────────────────────────────────────────────────────────

fn rate_matrix(n: usize) -> RateMatrix {
    RateMatrix::Dense(
        (0..n)
            .map(|i| (0..n).map(|j| GAMMA0 / (1.0 + D * (i as f64 - j as f64).abs())).collect())
            .collect(),
    )
}

/// Build initial state Σ_i Z_i for CoefficientThreshold.
fn initial_ct(n: usize, threshold: f64) -> PauliSum<Sct> {
    let mut p = PauliSum::builder().n_qubits(n).strategy(CoefficientThreshold(threshold)).build();
    let t = vec!['I'; n];
    for i in 0..n {
        let mut zi = t.clone();
        zi[i] = 'Z';
        p += (zi.iter().collect::<String>(), 1.0_f64);
    }
    p
}

/// Build initial state Σ_i Z_i for MaxPauliWeight.
fn initial_mpw(n: usize, max_weight: usize) -> PauliSum<Smpw> {
    let mut p = PauliSum::builder().n_qubits(n).strategy(MaxPauliWeight(max_weight)).build();
    let t = vec!['I'; n];
    for i in 0..n {
        let mut zi = t.clone();
        zi[i] = 'Z';
        p += (zi.iter().collect::<String>(), 1.0_f64);
    }
    p
}

/// Build initial state Σ_i Z_i for Budget.
fn initial_bud(n: usize, budget: Budget) -> PauliSum<Sbud> {
    let mut p = PauliSum::builder().n_qubits(n).strategy(budget).build();
    let t = vec!['I'; n];
    for i in 0..n {
        let mut zi = t.clone();
        zi[i] = 'Z';
        p += (zi.iter().collect::<String>(), 1.0_f64);
    }
    p
}

fn lindblad_ct(n: usize) -> LindbladOp<Sct> {
    let ops: Vec<JumpOp<Sct>> = (0..n)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
        .collect();
    LindbladOp::new(ops, rate_matrix(n))
}

fn lindblad_mpw(n: usize) -> LindbladOp<Smpw> {
    let ops: Vec<JumpOp<Smpw>> = (0..n)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
        .collect();
    LindbladOp::new(ops, rate_matrix(n))
}

fn lindblad_bud(n: usize) -> LindbladOp<Sbud> {
    let ops: Vec<JumpOp<Sbud>> = (0..n)
        .map(|i| JumpOp::Ladder(LadderOp { qubit: i, direction: LadderDirection::Lower }))
        .collect();
    LindbladOp::new(ops, rate_matrix(n))
}

/// Pattern for ⟨Z₀⟩: Z at position 0, I on all others.
/// Syntax: "Z0" = Z at index 0; the contains check rejects any word with
/// non-identity at other positions, giving the pure single-site Z₀ observable.
fn z0_pattern() -> PauliPattern {
    "Z0".into()
}

struct RunResult {
    wall_secs: f64,
    max_pauli: usize,
    obs:       f64,
}

fn run_ct(n: usize, threshold: f64, rtol: f64) -> RunResult {
    let lop  = lindblad_ct(n);
    let init = initial_ct(n, threshold);
    let pat  = z0_pattern();
    let t0 = Instant::now();
    let (_, out) = solve(
        None, &lop, &init,
        (0.0, TMAX), &[TMAX],
        |_, p: &PauliSum<Sct>| (p.trace(&pat), p.data().len()),
        SolverConfig { rtol, ..SolverConfig::default() },
    );
    RunResult {
        wall_secs: t0.elapsed().as_secs_f64(),
        max_pauli: out.iter().map(|&(_, sz)| sz).max().unwrap_or(0),
        obs:       out.last().map(|&(v, _)| v).unwrap_or(0.0),
    }
}

fn run_mpw(n: usize, max_weight: usize) -> RunResult {
    let lop  = lindblad_mpw(n);
    let init = initial_mpw(n, max_weight);
    let pat  = z0_pattern();
    let t0 = Instant::now();
    let (_, out) = solve(
        None, &lop, &init,
        (0.0, TMAX), &[TMAX],
        |_, p: &PauliSum<Smpw>| (p.trace(&pat), p.data().len()),
        SolverConfig::default(),
    );
    RunResult {
        wall_secs: t0.elapsed().as_secs_f64(),
        max_pauli: out.iter().map(|&(_, sz)| sz).max().unwrap_or(0),
        obs:       out.last().map(|&(v, _)| v).unwrap_or(0.0),
    }
}

fn run_bud(n: usize) -> RunResult {
    let budget = Budget { target: BUD_TARGET };
    let lop    = lindblad_bud(n);
    let init   = initial_bud(n, budget);
    let pat    = z0_pattern();
    let t0 = Instant::now();
    let (_, out) = solve(
        None, &lop, &init,
        (0.0, TMAX), &[TMAX],
        |_, p: &PauliSum<Sbud>| (p.trace(&pat), p.data().len()),
        // rtol matched to Budget min_threshold: absorbs the truncation-induced
        // k[6] perturbation and lets DOPRI5 take large steps.
        SolverConfig { rtol: BUD_RTOL, ..SolverConfig::default() },
    );
    RunResult {
        wall_secs: t0.elapsed().as_secs_f64(),
        max_pauli: out.iter().map(|&(_, sz)| sz).max().unwrap_or(0),
        obs:       out.last().map(|&(v, _)| v).unwrap_or(0.0),
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!(
        "n-atom superradiance scaling (T={TMAX}, γ₀={GAMMA0}, D={D})\n\
         Observable: ⟨Z₀⟩ at T={TMAX}\n"
    );
    println!(
        "{:<4} {:<22} {:>8} {:>12} {:>16}",
        "n", "Strategy", "|P|", "Time", "Infidelity"
    );
    println!("{}", "─".repeat(66));

    for &n in N_VALUES {
        // Reference run (treated as exact)
        // Reference uses default rtol (1e-6); tight threshold keeps state small and clean.
        let ref_result = run_ct(n, REF_THRESHOLD, SolverConfig::default().rtol);

        let infidelity = |obs: f64| (obs - ref_result.obs).abs();

        // Helper closure for a formatted row
        let row = |label: &str, r: &RunResult, is_ref: bool| {
            let inf_str = if is_ref {
                "—".to_string()
            } else {
                format!("{:.2e}", infidelity(r.obs))
            };
            println!(
                "{:<4} {:<22} {:>8} {:>12.3?} {:>16}",
                n, label, r.max_pauli,
                std::time::Duration::from_secs_f64(r.wall_secs),
                inf_str,
            );
        };

        row("Reference(1e-8)",  &ref_result, true);
        // CT_RTOL = 10×CT_THRESHOLD: prevents DOPRI5 fighting the pruning perturbation.
        row("CT(1e-3+rtol)",    &run_ct(n, CT_THRESHOLD, CT_RTOL), false);
        row("MPW(2)",           &run_mpw(n, MPW2_WEIGHT), false);
        row("MPW(4)",           &run_mpw(n, MPW4_WEIGHT), false);
        row(&format!("Budget({BUD_TARGET},1e-4)"), &run_bud(n), false);
        println!("{}", "─".repeat(66));
    }

    println!();
    println!("Notes:");
    println!("  Reference      : CoefficientThreshold({REF_THRESHOLD:.0e}), default rtol — treated as exact");
    println!("  CT(1e-3+rtol)  : CoefficientThreshold({CT_THRESHOLD:.0e}), rtol={CT_RTOL:.0e} (10×threshold, avoids DOPRI5 rejections)");
    println!("  MPW(w)         : MaxPauliWeight(w), default rtol — hard O(nʷ) bound; slow when n is large");
    println!("                   because DOPRI5 fights the weight-truncation k[6] perturbation (same cause");
    println!("                   as Budget without matched rtol); no natural threshold to derive rtol from");
    println!("  Budget(B,1e-4) : at most {BUD_TARGET} terms, rtol={BUD_RTOL:.0e} — always fast, accuracy degrades as cap bites");
    println!();
    println!("Polynomial scaling: MPW(w) keeps |P| ≤ Σ_{{k=0}}^{{w}} C(n,k)·3^k  (polynomial in n for fixed w).");
    println!("Budget:             |P| ≤ {BUD_TARGET} always — O(1) in n, but infidelity may grow as the cap becomes tight.");
}
