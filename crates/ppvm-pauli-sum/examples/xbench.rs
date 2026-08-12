// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The `ppvm` side of the cross-library Pauli-propagation benchmark.
//!
//! Two Trotter workloads — TFIM magnetization and a Heisenberg autocorrelator —
//! propagated in the Heisenberg picture with a coefficient-magnitude truncation
//! after every gate. See `benchmarks/cross-library/README.md` for the shared
//! spec: the gate order, the `θ = 2·c·dt` convention, the environment contract,
//! and the CSV schema that every runner in that harness prints.
//!
//! ```bash
//! MODEL=tfim QUBITS=8,16,24 STEPS=10 DT=0.1 JCOUP=1.0 HFIELD=1.0 ATOL=1e-6 \
//!   cargo run --release -p ppvm-pauli-sum --example xbench
//! ```

use std::time::Instant;

use ppvm_pauli_sum::prelude::*;
use ppvm_pauli_sum::strategy::CoefficientThreshold;

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

/// Smallest power-of-two byte width that holds `n` qubits with room to spare.
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

type Cfg<const N: usize> = config::fxhash::Byte<N, f64, CoefficientThreshold, PauliWord<[u8; N]>>;

/// Build the seed observable for `model` on `n` sites.
fn seed<const N: usize>(n: usize, p: Params) -> PauliSum<Cfg<N>> {
    let mut sum: PauliSum<Cfg<N>> = PauliSum::builder()
        .n_qubits(n)
        .strategy(CoefficientThreshold(p.atol))
        .capacity(1 << 12)
        .build();
    match p.model {
        Model::Tfim => {
            for i in 0..n {
                sum += (PauliWord::from(site_word(n, i, 'Z').as_str()), 1.0);
            }
        }
        Model::Heisenberg => {
            sum += (PauliWord::from(site_word(n, 0, 'Z').as_str()), 1.0);
        }
    }
    sum
}

/// Propagate the shared gate sequence through `state`.
fn propagate<const N: usize>(state: &mut PauliSum<Cfg<N>>, n: usize, p: Params) {
    let theta_bond = 2.0 * p.j * p.dt;
    let theta_site = 2.0 * p.h * p.dt;
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
}

/// `⟨0…0|O|0…0⟩` for TFIM; the `Z_0` coefficient for Heisenberg.
///
/// The diagonal contraction is spelled out rather than delegated: `⟨0|Z|0⟩ = 1`
/// and `⟨0|X|0⟩ = ⟨0|Y|0⟩ = 0`, so exactly the X-free terms survive and each
/// contributes its own coefficient.
fn readout<const N: usize>(state: &PauliSum<Cfg<N>>, n: usize, p: Params) -> f64 {
    match p.model {
        Model::Tfim => state
            .iter()
            .filter(|(word, _)| (0..n).all(|i| !word.get_xbit(i)))
            .map(|(_, c)| *c)
            .sum(),
        Model::Heisenberg => state
            .data()
            .get(&PauliWord::from(site_word(n, 0, 'Z').as_str()))
            .copied()
            .unwrap_or(0.0),
    }
}

/// Run one model at width `n`, returning `(best seconds, final support, observable)`.
fn run<const N: usize>(n: usize, p: Params) -> (f64, usize, f64) {
    let base = seed::<N>(n, p);
    let mut best = f64::INFINITY;
    let mut terms = 0usize;
    let mut observable = f64::NAN;
    for _ in 0..p.iters {
        let mut state = base.clone();
        let t0 = Instant::now();
        propagate(&mut state, n, p);
        best = best.min(t0.elapsed().as_secs_f64());
        terms = state.len();
        observable = readout(&state, n, p);
    }
    (best, terms, observable)
}

/// Print the whole propagated support as `word coefficient`, largest first —
/// the format the driver diffs across engines.
fn dump<const N: usize>(n: usize, p: Params) {
    let mut state = seed::<N>(n, p);
    propagate(&mut state, n, p);
    let mut out: Vec<(String, f64)> = state.iter().map(|(k, c)| (k.to_string(), *c)).collect();
    out.sort_by(|a, b| {
        b.1.abs()
            .partial_cmp(&a.1.abs())
            .unwrap()
            .then(a.0.cmp(&b.0))
    });
    println!("# {} terms", out.len());
    for (word, coeff) in out {
        println!("{word} {coeff:+.12e}");
    }
}

macro_rules! dispatch {
    ($f:ident, $n:expr, $p:expr) => {
        match storage_bytes($n) {
            2 => $f::<2>($n, $p),
            4 => $f::<4>($n, $p),
            8 => $f::<8>($n, $p),
            16 => $f::<16>($n, $p),
            32 => $f::<32>($n, $p),
            64 => $f::<64>($n, $p),
            b => panic!("no storage tier for {b} bytes (n = {})", $n),
        }
    };
}

fn main() {
    let p = Params::from_env();
    let qubits: Vec<usize> = std::env::var("QUBITS")
        .ok()
        .map(|s| s.split(',').filter_map(|t| t.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![8, 12, 16, 20, 24, 28, 32]);

    if std::env::var("DUMP").is_ok() {
        let n = qubits[0];
        dispatch!(dump, n, p);
        return;
    }

    eprintln!(
        "ppvm {}: steps={} dt={} J={} h={} atol={:e} iters={}",
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
        let (t, terms, obs) = dispatch!(run, n, p);
        println!(
            "{},ppvm,{},{},{},{:e},{:.6},{},{:.12e}",
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
