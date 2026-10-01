// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared CLI, timing and JSON output for `trotter_cpu` / `trotter_gpu`.
//!
//! Defaults reproduce `ppvm-benchmarks/trotter-benchmarks`; the JSON follows
//! that repo's `results*.json` schema so `plot.py` can read it.

use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use clap::{Args, ValueEnum};
use ppvm_cupauliprop::Propagator;
use ppvm_cupauliprop::circuits::{heisenberg, lattice_edges, sum_z_terms, tfim};
use serde_json::{Value, json};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Model {
    /// 1D transverse-field Ising chain (`test_trotter.py`).
    Tfim,
    /// Heisenberg model on an L×L open-BC lattice (`test_heisenberg.py`).
    Heisenberg,
}

#[derive(Debug, Args)]
pub struct BenchArgs {
    #[arg(long, value_enum, default_value = "tfim")]
    pub model: Model,
    /// System sizes: N for tfim, L for heisenberg [default: the ppvm-benchmarks scan].
    #[arg(long, value_delimiter = ',')]
    pub sizes: Option<Vec<usize>>,
    /// Timed rounds per size [default: 10 for tfim, 3 for heisenberg].
    #[arg(long)]
    pub rounds: Option<usize>,
    #[arg(long, default_value_t = 20)]
    pub steps: usize,
    #[arg(long, default_value_t = 0.1)]
    pub dt: f64,
    /// Coupling J.
    #[arg(long, default_value_t = 1.0)]
    pub j: f64,
    /// Transverse field h (tfim only).
    #[arg(long, default_value_t = 1.0)]
    pub h: f64,
    /// Depolarizing parameter p; applied as `pauli_error(q, [p/4; 3])`.
    #[arg(long, default_value_t = 1e-4)]
    pub noise: f64,
    /// Coefficient truncation threshold (`min_abs_coeff`).
    #[arg(long, default_value_t = 1e-6)]
    pub cutoff: f64,
    /// Write results JSON here.
    #[arg(long)]
    pub out: Option<PathBuf>,
}

impl BenchArgs {
    pub fn sizes(&self) -> Vec<usize> {
        self.sizes.clone().unwrap_or_else(|| match self.model {
            Model::Tfim => vec![4, 8, 12, 16, 24, 32, 48, 64, 80, 100, 128],
            Model::Heisenberg => vec![2, 3, 4, 5, 6],
        })
    }

    pub fn rounds(&self) -> usize {
        self.rounds.unwrap_or(match self.model {
            Model::Tfim => 10,
            Model::Heisenberg => 3,
        })
    }

    pub fn n_qubits(&self, size: usize) -> usize {
        match self.model {
            Model::Tfim => size,
            Model::Heisenberg => size * size,
        }
    }

    fn run_circuit<P: Propagator>(&self, state: &mut P, size: usize) {
        let noise = [self.noise / 4.0; 3];
        match self.model {
            Model::Tfim => tfim(
                state,
                size,
                self.steps,
                self.dt * self.h,
                self.dt * self.j,
                noise,
            ),
            Model::Heisenberg => heisenberg(
                state,
                &lattice_edges(size),
                self.steps,
                self.dt * self.j,
                noise,
            ),
        }
    }
}

/// Time `rounds` propagations of one size. `build` creates a fresh state
/// (untimed); `finish` runs inside the timed region after the circuit
/// (e.g. a device sync).
pub fn bench_size<P: Propagator>(
    args: &BenchArgs,
    size: usize,
    build: impl Fn(usize, &[String]) -> P,
    finish: impl Fn(&P),
) -> Value {
    let n_qubits = args.n_qubits(size);
    let terms = sum_z_terms(n_qubits);
    let mut times = Vec::with_capacity(args.rounds());
    let mut last = None;
    for _ in 0..args.rounds() {
        // Free the previous state first: GPU states can reserve most of device memory.
        drop(last.take());
        let mut state = build(n_qubits, &terms);
        let t0 = Instant::now();
        args.run_circuit(&mut state, size);
        finish(&state);
        times.push(t0.elapsed().as_secs_f64());
        last = Some(state);
    }
    let state = last.expect("at least one round");
    let expectation = state.overlap_with_zero();
    let n_terms = state.len();

    let mut sorted = times.clone();
    sorted.sort_by(f64::total_cmp);
    let min = sorted[0];
    let mean = times.iter().sum::<f64>() / times.len() as f64;
    let median = if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2]
    } else {
        0.5 * (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2])
    };
    eprintln!(
        "n_qubits={n_qubits:>4}  min={min:.6}s  median={median:.6}s  rounds={}  \
         <O>={expectation:.15}  terms={n_terms}",
        times.len()
    );

    let mut row = json!({
        "n_qubits": n_qubits,
        "n_samples": times.len(),
        "time_seconds_min": min,
        "time_seconds_mean": mean,
        "time_seconds_median": median,
        "expectation_value": expectation,
        "n_terms_final": n_terms,
    });
    if args.model == Model::Heisenberg {
        row["L"] = json!(size);
    }
    row
}

fn git_commit() -> Option<String> {
    let out = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_owned())
}

/// Assemble the results document and write it to `--out` (if given).
pub fn write_results(args: &BenchArgs, library: &str, extra_info: Value, results: Vec<Value>) {
    let mut params = json!({
        "dt": args.dt,
        "time": args.dt * args.steps as f64,
        "n_steps": args.steps,
        "j": args.j,
        "noise_param": args.noise,
        "min_abs_coeff": args.cutoff,
        "observable": "sum_i Z_i",
    });
    match args.model {
        Model::Tfim => params["h"] = json!(args.h),
        Model::Heisenberg => {
            params["lattice"] = json!("2D square, open BC");
            params["hamiltonian"] = json!("Heisenberg: J Σ_{<i,j>} (X_iX_j + Y_iY_j + Z_iZ_j)");
        }
    }
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let doc = json!({
        "library": library,
        "run_info": {
            "timestamp_unix": timestamp,
            "git_commit": git_commit(),
            "command": std::env::args().collect::<Vec<_>>().join(" "),
            "rayon_num_threads": std::env::var("RAYON_NUM_THREADS").ok(),
            "rounds": args.rounds(),
            "backend": extra_info,
        },
        "params": params,
        "results": results,
    });
    if let Some(path) = &args.out {
        let text = serde_json::to_string_pretty(&doc).expect("serialize results");
        std::fs::write(path, text + "\n").expect("write results file");
        eprintln!("wrote {}", path.display());
    }
}
