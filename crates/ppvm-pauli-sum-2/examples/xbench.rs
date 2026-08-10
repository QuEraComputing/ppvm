// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The `ppvm-*-2` side of the cross-library Pauli-propagation benchmark.
//!
//! Two models, both driven as an explicit first-order Trotter product of Pauli
//! rotations so that every library in the comparison executes the *same gate
//! list in the same order* rather than its own idea of a Trotter step:
//!
//! * `tfim` — `H = J Σ Z_i Z_{i+1} + h Σ X_i`. Per step: `RX(2h·dt)` on every
//!   site in ascending order, then `RZZ(2J·dt)` on every bond in ascending
//!   order. Observable `O = Σ_i Z_i`, read out as `⟨0…0|O(t)|0…0⟩`.
//! * `heisenberg` — `H = J Σ (X_iX_{i+1} + Y_iY_{i+1} + Z_iZ_{i+1}) + h Σ Z_i`.
//!   Per step: `RXX`, `RYY`, `RZZ` at `2J·dt` on every bond in ascending order,
//!   then `RZ(2h·dt)` on every site. Observable `O = Z_0`, read out as the
//!   autocorrelator `S(t) = tr[Z_0·O(t)]/2^n`, i.e. the coefficient of `Z_0`.
//!
//! `θ = 2·c·dt` for a Hamiltonian term `c·G` is the convention every engine in
//! the comparison uses for `exp(iθ/2·G)·P·exp(−iθ/2·G)`, so the propagated
//! operator is identical up to floating-point associativity — the driver
//! asserts that.
//!
//! Truncation is a coefficient-magnitude threshold applied **after every gate**,
//! which is the one truncation rule all four libraries share.
//!
//! Usage — parameters come from the environment so every runner reads the same
//! contract; CSV goes to stdout and progress to stderr:
//!
//! ```bash
//! MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
//!   cargo run --release -p ppvm-pauli-sum-2 --example xbench
//! ```

use std::time::Instant;

use ppvm_pauli_sum_2::{CoefficientThreshold, HashMapStore, PauliPattern, PauliWord, Sum};
use ppvm_traits_2::{RotationOne, RotationTwo, Trace};

/// The model, its Hamiltonian couplings, and the Trotter schedule.
#[derive(Clone, Copy)]
struct Params {
    model: Model,
    steps: usize,
    dt: f64,
    j: f64,
    h: f64,
    atol: f64,
    iters: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Model {
    Tfim,
    Heisenberg,
}

impl Model {
    fn parse(s: &str) -> Self {
        match s {
            "tfim" => Model::Tfim,
            "heisenberg" => Model::Heisenberg,
            other => panic!("unknown MODEL {other:?} (expected `tfim` or `heisenberg`)"),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Model::Tfim => "tfim",
            Model::Heisenberg => "heisenberg",
        }
    }
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

impl Params {
    fn from_env() -> Self {
        Params {
            model: Model::parse(&std::env::var("MODEL").unwrap_or_else(|_| "tfim".to_string())),
            steps: env_usize("STEPS", 10),
            dt: env_f64("DT", 0.1),
            j: env_f64("JCOUP", 1.0),
            h: env_f64("HFIELD", 1.0),
            atol: env_f64("ATOL", 1e-6),
            iters: env_usize("ITERS", 3),
        }
    }
}

/// Smallest power-of-two byte width that holds `n` qubits with room to spare —
/// the same storage-tier ladder the Python bindings dispatch on.
fn storage_bytes(n: usize) -> usize {
    let need = n.div_ceil(8);
    let mut k = 0;
    while (1usize << k) <= need {
        k += 1;
    }
    1usize << k
}

/// A single-site Pauli word `p` at `site`, identity elsewhere.
fn site_word(n: usize, site: usize, p: char) -> String {
    (0..n).map(|j| if j == site { p } else { 'I' }).collect()
}

/// Run one model at width `n`, returning `(best seconds, final support, observable)`.
fn run<const B: usize>(n: usize, p: Params) -> (f64, usize, f64) {
    type Store<const B: usize> = HashMapStore<PauliWord<[u8; B]>, f64>;
    type SumOf<const B: usize> = Sum<Store<B>, CoefficientThreshold>;

    let policy = CoefficientThreshold { threshold: p.atol };
    // The seed observable, built once and cloned per timed iteration.
    let mut seed: SumOf<B> = Sum::with_capacity(n, policy, 1 << 12);
    match p.model {
        // Total magnetization Σ_i Z_i.
        Model::Tfim => {
            for i in 0..n {
                seed += (PauliWord::from(site_word(n, i, 'Z').as_str()), 1.0);
            }
        }
        // A single site, whose autocorrelator is the readout.
        Model::Heisenberg => {
            seed += (PauliWord::from(site_word(n, 0, 'Z').as_str()), 1.0);
        }
    }

    let theta_bond = 2.0 * p.j * p.dt;
    let theta_site = 2.0 * p.h * p.dt;

    let mut best = f64::INFINITY;
    let mut terms = 0usize;
    let mut observable = f64::NAN;
    for _ in 0..p.iters {
        let mut state = seed.clone();
        let t0 = Instant::now();
        for _ in 0..p.steps {
            match p.model {
                Model::Tfim => {
                    for i in 0..n {
                        state.rx(i, theta_site);
                        state.truncate();
                    }
                    for i in 0..n.saturating_sub(1) {
                        state.rzz(i, i + 1, theta_bond);
                        state.truncate();
                    }
                }
                Model::Heisenberg => {
                    for i in 0..n.saturating_sub(1) {
                        state.rxx(i, i + 1, theta_bond);
                        state.truncate();
                        state.ryy(i, i + 1, theta_bond);
                        state.truncate();
                        state.rzz(i, i + 1, theta_bond);
                        state.truncate();
                    }
                    for i in 0..n {
                        state.rz(i, theta_site);
                        state.truncate();
                    }
                }
            }
        }
        best = best.min(t0.elapsed().as_secs_f64());
        terms = state.len();
        observable = match p.model {
            // ⟨0…0|O|0…0⟩ — the `Z?*` contraction keeps exactly the terms with
            // no X/Y factor, each of which evaluates to its own coefficient.
            Model::Tfim => Trace::trace(&state, &PauliPattern::zero_state()),
            // tr[Z_0·O(t)]/2^n is the coefficient of Z_0, the Paulis being
            // orthonormal under that pairing.
            Model::Heisenberg => state
                .get(&PauliWord::from(site_word(n, 0, 'Z').as_str()))
                .unwrap_or(0.0),
        };
    }
    (best, terms, observable)
}

/// Print the whole propagated support as `word coefficient`, largest first.
///
/// The cross-library driver diffs this against the other engines' dumps: an
/// observable is one number and can agree by luck, but a term-for-term match
/// says the four engines really did propagate the same operator. Small widths
/// only — this is a debugging/validation path, not part of the timed run.
fn dump<const B: usize>(n: usize, p: Params) {
    let policy = CoefficientThreshold { threshold: p.atol };
    let mut state: Sum<HashMapStore<PauliWord<[u8; B]>, f64>, CoefficientThreshold> =
        Sum::with_capacity(n, policy, 1 << 12);
    match p.model {
        Model::Tfim => {
            for i in 0..n {
                state += (PauliWord::from(site_word(n, i, 'Z').as_str()), 1.0);
            }
        }
        Model::Heisenberg => {
            state += (PauliWord::from(site_word(n, 0, 'Z').as_str()), 1.0);
        }
    }
    let (theta_bond, theta_site) = (2.0 * p.j * p.dt, 2.0 * p.h * p.dt);
    for _ in 0..p.steps {
        match p.model {
            Model::Tfim => {
                for i in 0..n {
                    state.rx(i, theta_site);
                    state.truncate();
                }
                for i in 0..n.saturating_sub(1) {
                    state.rzz(i, i + 1, theta_bond);
                    state.truncate();
                }
            }
            Model::Heisenberg => {
                for i in 0..n.saturating_sub(1) {
                    state.rxx(i, i + 1, theta_bond);
                    state.truncate();
                    state.ryy(i, i + 1, theta_bond);
                    state.truncate();
                    state.rzz(i, i + 1, theta_bond);
                    state.truncate();
                }
                for i in 0..n {
                    state.rz(i, theta_site);
                    state.truncate();
                }
            }
        }
    }
    let mut terms: Vec<(String, f64)> = Vec::with_capacity(state.len());
    state.for_each_ref(|k, c| terms.push((k.to_string(), *c)));
    terms.sort_by(|a, b| {
        b.1.abs()
            .partial_cmp(&a.1.abs())
            .unwrap()
            .then(a.0.cmp(&b.0))
    });
    println!("# {} terms", terms.len());
    for (word, coeff) in terms {
        println!("{word} {coeff:+.12e}");
    }
}

fn main() {
    let p = Params::from_env();

    // Validation path: dump the propagated support instead of timing it.
    if std::env::var("DUMP").is_ok() {
        let n: usize = std::env::var("QUBITS")
            .ok()
            .and_then(|s| s.split(',').next().and_then(|t| t.trim().parse().ok()))
            .unwrap_or(4);
        match storage_bytes(n) {
            2 => dump::<2>(n, p),
            4 => dump::<4>(n, p),
            8 => dump::<8>(n, p),
            16 => dump::<16>(n, p),
            32 => dump::<32>(n, p),
            64 => dump::<64>(n, p),
            b => panic!("no storage tier for {b} bytes (n = {n})"),
        }
        return;
    }

    let qubits: Vec<usize> = std::env::var("QUBITS")
        .ok()
        .map(|s| s.split(',').filter_map(|t| t.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![8, 12, 16, 20, 24, 28, 32]);

    eprintln!(
        "ppvm-2 {}: steps={} dt={} J={} h={} atol={:e} iters={}",
        p.model.name(),
        p.steps,
        p.dt,
        p.j,
        p.h,
        p.atol,
        p.iters
    );

    println!("model,library,qubits,steps,dt,atol,time_s,terms,observable");
    for &n in &qubits {
        let (t, terms, obs) = match storage_bytes(n) {
            2 => run::<2>(n, p),
            4 => run::<4>(n, p),
            8 => run::<8>(n, p),
            16 => run::<16>(n, p),
            32 => run::<32>(n, p),
            64 => run::<64>(n, p),
            b => panic!("no storage tier for {b} bytes (n = {n})"),
        };
        println!(
            "{},ppvm-2,{},{},{},{:e},{:.6},{},{:.12e}",
            p.model.name(),
            n,
            p.steps,
            p.dt,
            p.atol,
            t,
            terms,
            obs
        );
        eprintln!("  n={n:3}  {t:9.4}s  {terms:>9} terms  obs={obs:+.9e}");
        use std::io::Write;
        std::io::stdout().flush().ok();
    }
}
