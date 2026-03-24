use std::array;

use ppvm_runtime::{config::fxhash::ByteF64, prelude::*, strategy::CoefficientThreshold};
use ppvm_timeevolve::{Budget, CollapseOp, LindbladOp, RateMatrix, SolverConfig, solve::solve};

const NBYTES: usize = 2;
type S = ByteF64<NBYTES, CoefficientThreshold>;
// Budget variant: caps |P| at 150 entries while matching rtol to the truncation threshold.
// When |P| ≤ 150, behaviour is identical to CoefficientThreshold(BUDGET_THRESHOLD).
type SB = ByteF64<NBYTES, Budget>;
const BUDGET_THRESHOLD: f64 = 1e-6;
const BUDGET_TARGET: usize = 150;

fn build_lindblad_ops<T: ppvm_runtime::config::Config>(
    n: usize,
    gamma_mat: RateMatrix,
    phase_y: u8,
) -> LindbladOp<T>
where
    ppvm_runtime::prelude::PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>:
        std::ops::Mul<Output = ppvm_runtime::prelude::PhasedPauliWord<T::Storage, T::BuildHasher, T::PauliWordType>>
        + std::ops::MulAssign
        + Clone,
    T::PauliWordType: Clone
        + for<'a> From<&'a str>
        + std::borrow::Borrow<ppvm_runtime::prelude::PauliWord<T::Storage, T::BuildHasher>>,
{
    let ppw = |pauli: &str, phase: u8| {
        ppvm_runtime::prelude::PhasedPauliWord::<T::Storage, T::BuildHasher, T::PauliWordType>
            ::build_from_word(T::PauliWordType::from(pauli), phase)
    };
    let tmp = vec!['I'; n];
    let mut c_ops: Vec<CollapseOp<T>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut c = CollapseOp::<T>::new(n);
        let mut px = tmp.clone();
        let mut py = tmp.clone();
        px[i] = 'X';
        py[i] = 'Y';
        c.push(ppw(&px.iter().collect::<String>(), 0), 1.0);
        c.push(ppw(&py.iter().collect::<String>(), phase_y), 1.0);
        c_ops.push(c);
    }
    LindbladOp::new(c_ops, gamma_mat)
}

fn main() {
    let n = 5;
    let gamma0 = 1.0;
    let d = 0.1;
    let tmax = 1.0;
    const TSTEPS: usize = 100;
    let dt = tmax / TSTEPS as f64;

    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(n);
    for i in 0i32..(n as i32) {
        let mut row: Vec<f64> = Vec::with_capacity(n);
        for j in 0i32..(n as i32) {
            row.push(if i == j { gamma0 } else { gamma0 / (1.0 + d * (i - j).unsigned_abs() as f64) });
        }
        rows.push(row);
    }

    let zero_pattern: PauliPattern = "Z?*".into();
    let save_at: [f64; TSTEPS] = array::from_fn(|i| dt * i as f64);
    let tmp = vec!['I'; n];

    // ── Baseline: CoefficientThreshold(1e-6) ────────────────────────────────
    let t0 = std::time::Instant::now();
    let gamma_mat_base = RateMatrix::Dense(rows.clone());
    let lindblad_base = build_lindblad_ops::<S>(n, gamma_mat_base, 3);
    let strat_base = CoefficientThreshold(1e-6);
    let mut initial_base: PauliSum<S> =
        PauliSum::builder().n_qubits(n).strategy(strat_base).build();
    for i in 0..n {
        let mut zi = tmp.clone();
        zi[i] = 'Z';
        initial_base += (zi.iter().collect::<String>(), 1.0);
    }
    let (_, baseline_vals) = solve(
        None, &lindblad_base, &initial_base, (0.0, tmax), &save_at,
        |_, p: &PauliSum<S>| p.trace(&zero_pattern),
        SolverConfig::default(),
    );
    let elapsed_base = t0.elapsed();

    // ── Budget: Budget { target=150, min_threshold=1e-6 } + matched rtol ────
    let t1 = std::time::Instant::now();
    let gamma_mat_bud = RateMatrix::Dense(rows.clone());
    let lindblad_bud = build_lindblad_ops::<SB>(n, gamma_mat_bud, 3);
    let strat_bud = Budget { target: BUDGET_TARGET, min_threshold: BUDGET_THRESHOLD };
    let mut initial_bud: PauliSum<SB> =
        PauliSum::builder().n_qubits(n).strategy(strat_bud).build();
    for i in 0..n {
        let mut zi = tmp.clone();
        zi[i] = 'Z';
        initial_bud += (zi.iter().collect::<String>(), 1.0);
    }
    let (_, budget_vals) = solve(
        None, &lindblad_bud, &initial_bud, (0.0, tmax), &save_at,
        |_, p: &PauliSum<SB>| p.trace(&zero_pattern),
        // Match rtol to truncation threshold — DOPRI5 takes larger steps when smooth.
        SolverConfig { rtol: BUDGET_THRESHOLD, ..SolverConfig::default() },
    );
    let elapsed_bud = t1.elapsed();

    // ── Report ───────────────────────────────────────────────────────────────
    println!("Gamma: {:?}", rows);
    println!();
    println!("{:<12} {:>12} {:>10} {:>14}", "Variant", "Wall time", "Speedup", "Max |fidelity err|");
    println!("{}", "-".repeat(52));

    let fidelity_err: f64 = baseline_vals.iter().zip(budget_vals.iter())
        .map(|(b, u)| (b - u).abs())
        .fold(0.0_f64, f64::max);
    let speedup = elapsed_base.as_secs_f64() / elapsed_bud.as_secs_f64();

    println!("{:<12} {:>12.3?} {:>10} {:>14}", "Baseline", elapsed_base, "1.00×", "—");
    println!("{:<12} {:>12.3?} {:>9.2}× {:>14.2e}", "Budget", elapsed_bud, speedup, fidelity_err);
    println!();
    println!("tout: {:?}", save_at);
    println!("baseline values: {:?}", baseline_vals);
    println!("budget   values: {:?}", budget_vals);
}
